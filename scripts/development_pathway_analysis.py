"""Standalone development-candidate analysis; writes only its dedicated directories."""
from pathlib import Path
from itertools import combinations
import hashlib
import importlib.metadata
import json
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from rnaseq_inputs import ROOT, write_json, sha256, read_counts, read_metadata, align_metadata, read_contrast
from development_pathway_matching import load_rules, read_annotations, match_annotations, build_memberships, target_counts, naming_status

OUT=ROOT/'results/processedData/development_pathways'
FIG=ROOT/'results/figures/development_pathways'
NOTEBOOK=ROOT/'scripts/05_development_pathway_analysis.ipynb'
NEW_SCRIPTS={'development_pathway_analysis.py','development_pathway_matching.py','development_pathway_smoke_tests.py',
             'development_pathway_config.json','05_development_pathway_analysis.ipynb'}


def save_table(table,path,index=False):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.name.endswith('.gz'):
        table.to_csv(path,sep='\t',index=index,float_format='%.17g',compression={'method':'gzip','mtime':0})
    else:
        table.to_csv(path,sep='\t',index=index,float_format='%.17g')


def settings(overrides=None):
    cfg=json.loads((ROOT/'scripts/development_pathway_config.json').read_text())
    cfg.update(overrides or {})
    if not 0<cfg['padj_threshold']<=1 or cfg['lfc_threshold']<0:
        raise ValueError('Invalid DE thresholds')
    return cfg


def slug(target):
    return re.sub(r'[^A-Za-z0-9_.-]+','_',target).strip('_')


def savefig(fig,name):
    FIG.mkdir(parents=True,exist_ok=True)
    for extension in ['png','pdf']:
        metadata={'CreationDate':None,'ModDate':None} if extension=='pdf' else {'Software':'HydraAEP development candidate analysis'}
        fig.savefig(FIG/f'{name}.{extension}',dpi=170,bbox_inches='tight',metadata=metadata)
    plt.close(fig)


def protected_paths():
    paths=[p for p in (ROOT/'data').rglob('*') if p.is_file()]
    paths += [p for p in (ROOT/'results').rglob('*') if p.is_file() and not p.is_relative_to(OUT) and not p.is_relative_to(FIG)]
    paths += [p for p in (ROOT/'scripts').iterdir() if p.is_file() and p.name not in NEW_SCRIPTS]
    paths.append(ROOT/'requirements.txt')
    return sorted(set(paths))


def start_analysis(overrides=None):
    cfg=settings(overrides)
    OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    manifest=OUT/'preserved_files_sha256.tsv'
    if manifest.exists():
        previous=pd.read_csv(manifest,sep='\t')
        changed=[r.path for r in previous.itertuples() if not (ROOT/r.path).is_file() or sha256(ROOT/r.path)!=r.sha256]
        if changed:
            raise ValueError(f'Protected inputs/existing outputs changed: {changed}')
    else:
        save_table(pd.DataFrame([dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha256(p)) for p in protected_paths()]),manifest)
    write_json(OUT/'configuration_used.json',cfg)
    requested=ROOT/'results/deseq2_rsem/tables'
    external=json.loads((ROOT/cfg['contrast_manifest_file']).read_text())
    items=external['contrasts']
    if len({item['name'] for item in items})!=len(items):
        raise ValueError('Duplicate DE source labels')
    tables={};issues=[]
    for item in items:
        try:
            tables[item['name']]=read_contrast(item)
        except (ValueError,OSError,KeyError) as error:
            issues.append(dict(source=item.get('name'),error=str(error)))
    if not tables:
        issues.append(dict(step='DE',error='No validated existing DESeq2 source table'))
    raw,audit=read_annotations(ROOT/cfg['annotation_file'])
    rules=load_rules(ROOT/cfg['candidate_file'])
    if rules.target_group.map(slug).nunique()!=rules.target_group.nunique():
        raise ValueError('Target filename slugs collide')
    counts=files=meta=mapping=None
    try:
        meta,controls=read_metadata()
        counts,files,transcripts=read_counts()
        mapping=align_metadata(files,meta,'')
        save_table(files,OUT/'count_file_audit.tsv')
        save_table(mapping,OUT/'sample_file_mapping.tsv')
        save_table(controls,OUT/'control_metadata.tsv')
        save_table(meta.groupby(['EXP','treatment','batch']).size().rename('metadata_records').reset_index(),OUT/'metadata_design_counts.tsv')
    except (ValueError,OSError,KeyError) as error:
        issues.append(dict(step='sample_expression',error=str(error)))
    original_DE_headers={item['name']:list(pd.read_csv(ROOT/item['path'],sep=item.get('delimiter'),engine='python',nrows=0)) for item in items if item['name'] in tables}
    headers=dict(candidate=list(pd.read_csv(ROOT/cfg['candidate_file'],nrows=0)),annotation=audit['header'],DE=original_DE_headers,
                 derived_candidate_columns=list(rules),derived_DE_columns={name:list(table.reset_index()) for name,table in tables.items()},
                 metadata=list(meta) if meta is not None else [],count_headers=sorted(files['columns'].unique()) if files is not None else [],
                 requested_DE_directory=requested.relative_to(ROOT).as_posix(),requested_DE_directory_exists=requested.exists(),
                 verified_DE_paths=[item['path'] for item in items if item['name'] in tables])
    write_json(OUT/'input_headers_and_paths.json',headers)
    write_json(OUT/'contrasts_used.json',items)
    write_json(OUT/'annotation_parser_audit.json',audit)
    write_json(OUT/'input_issues.json',issues)
    save_table(raw.loc[raw['query'].duplicated(keep=False)],OUT/'duplicated_annotation_queries.tsv')
    status=dict(sample_treatment_expression='blocked',reason='RSEM contents have no internal sample identifier; exact filename/metadata ID pairs are candidates only.',
                file_expression='available' if counts is not None else 'blocked',
                filename_candidate_matches=int(mapping.mapping_status.eq('candidate_filename_match').sum()) if mapping is not None else None,
                files_without_metadata=int(mapping.mapping_status.eq('file_without_metadata').sum()) if mapping is not None else None,
                metadata_without_files=int(mapping.mapping_status.eq('metadata_without_file').sum()) if mapping is not None else None,
                biological_DE_directions='unresolved for sources whose direction_verified flag is false; up/down lists refer only to supplied LFC signs')
    write_json(OUT/'expression_mapping_status.json',status)
    print(f"Inspected {len(rules)} rules, {rules.target_group.nunique()} targets, {len(raw)} annotation records, {len(tables)} DE tables",flush=True)
    return dict(cfg=cfg,rules=rules,raw=raw,annotation_audit=audit,tables=tables,counts=counts,files=files,meta=meta,mapping=mapping,issues=issues,expression_status=status)


def extract_candidates(state):
    cfg,rules,raw=state['cfg'],state['rules'],state['raw']
    evidence,expansions=match_annotations(rules,raw,cfg['supplementary_OG_matching'])
    membership,catalog=build_memberships(evidence,raw)
    summary=target_counts(rules,evidence,membership)
    save_table(rules,OUT/'candidate_rules.tsv')
    save_table(expansions,OUT/'synonym_expansions_and_patterns.tsv')
    save_table(evidence,OUT/'gene_match_evidence.tsv.gz')
    save_table(membership,OUT/'gene_target_membership.tsv')
    save_table(summary,OUT/'target_gene_counts.tsv')
    save_table(catalog.loc[catalog['query'].isin(membership['query'])],OUT/'matched_gene_catalog.tsv')
    by_type=[]
    for target in rules.target_group.drop_duplicates():
        for kind in ['GO','KEGG_reference_pathway','gene_family_keyword']:
            subset=evidence.loc[evidence.target_group.eq(target) & evidence.annotation_type.eq(kind)]
            by_type.append(dict(target_group=target,annotation_type=kind,genes=subset['query'].nunique(),evidence_rows=len(subset),
                                primary_genes=subset.loc[subset.matched_field.ne('eggNOG_OGs'),'query'].nunique(),supplementary_OG_genes=subset.loc[subset.matched_field.eq('eggNOG_OGs'),'query'].nunique()))
    save_table(pd.DataFrame(by_type),OUT/'target_evidence_type_counts.tsv')
    per_rule=evidence.groupby('rule_id')['query'].nunique().rename('matched_genes')
    rule_counts=rules.merge(per_rule,on='rule_id',how='left',validate='one_to_one')
    rule_counts['matched_genes']=rule_counts.matched_genes.fillna(0).astype(int)
    save_table(rule_counts,OUT/'rule_match_counts.tsv')
    ann_ids=set(raw['query']);count_ids=set(state['counts'].index) if state['counts'] is not None else set()
    diagnostic=[];unmatched=[]
    sets={'annotations':ann_ids,'matched_candidates':set(membership['query'])}
    if state['counts'] is not None:
        sets['counts']=count_ids
    sets.update({f'DE:{name}':set(table.index) for name,table in state['tables'].items()})
    for source,ids in sets.items():
        diagnostic.append(dict(source=source,unique_ids=len(ids),duplicate_ids=0,matching_annotations=len(ids & ann_ids),
                               absent_from_annotations=len(ids-ann_ids),matching_counts=len(ids & count_ids) if count_ids else None,
                               absent_from_counts=len(ids-count_ids) if count_ids else None,mapping='exact HVAEP<scaffold>.G IDs; no suffix rewriting or positional joins'))
        for query in sorted(ids-ann_ids):
            unmatched.append(dict(source=source,query=query,reason='no_annotation_record'))
        if count_ids:
            for query in sorted(ids-count_ids):
                unmatched.append(dict(source=source,query=query,reason='absent_from_count_rows'))
    save_table(pd.DataFrame(diagnostic),OUT/'identifier_mapping_diagnostics.tsv')
    save_table(pd.DataFrame(unmatched,columns=['source','query','reason']),OUT/'unmatched_identifiers.tsv')
    gene_sets={target:set(membership.loc[membership.target_group.eq(target),'query']) for target in rules.target_group.drop_duplicates()}
    overlap=[]
    for a,b in combinations(gene_sets,2):
        shared=gene_sets[a]&gene_sets[b];union=gene_sets[a]|gene_sets[b]
        overlap.append(dict(target_a=a,target_b=b,shared_genes=len(shared),jaccard=len(shared)/len(union) if union else np.nan,genes=';'.join(sorted(shared))))
    save_table(pd.DataFrame(overlap),OUT/'target_overlap_pairs.tsv')
    manifest=[]
    for target,genes in gene_sets.items():
        target_slug=slug(target)
        directory=OUT/'gene_lists';directory.mkdir(exist_ok=True)
        path=directory/f'{target_slug}.txt'
        path.write_text(''.join(gene+'\n' for gene in sorted(genes)),encoding='utf-8')
        save_table(membership.loc[membership.target_group.eq(target)],directory/f'{target_slug}.tsv')
        manifest.append(dict(path=path.relative_to(OUT).as_posix(),target_group=target,contrast='',role='all_candidates',genes=len(genes),direction_basis='not applicable',direction_verified=False))
    state.update(evidence=evidence,expansions=expansions,membership=membership,catalog=catalog,summary=summary,gene_sets=gene_sets,list_manifest=manifest)
    print(f'Matched {membership["query"].nunique()} distinct genes in {len(membership)} gene/target memberships; {len(evidence)} evidence rows',flush=True)
    return summary


def analyze_de(state):
    cfg,membership=state['cfg'],state['membership']
    all_results=[];summaries=[];list_manifest=state['list_manifest']
    for name,table in state['tables'].items():
        de=table.reset_index().rename(columns={'gene_id':'query'})
        rows=membership.merge(de,on='query',how='left',validate='many_to_one',indicator=True)
        rows['contrast']=name
        rows['in_DE_export']=rows['_merge'].eq('both');rows=rows.drop(columns='_merge')
        rows['tested']=rows.tested.fillna(False).astype(bool)
        rows['eligible']=rows.eligible.fillna(False).astype(bool)
        rows['significant']=rows.eligible & rows.padj.lt(cfg['padj_threshold']) & rows.log2FoldChange.abs().ge(cfg['lfc_threshold'])
        rows['source_direction']=np.select([rows.significant & rows.log2FoldChange.gt(0),rows.significant & rows.log2FoldChange.lt(0)],['positive_LFC','negative_LFC'],default='not_significant_or_NA')
        rows['DE_status']=np.select([~rows.in_DE_export,~rows.tested,rows.padj.isna(),rows.log2FoldChange.isna(),rows.significant],
                                    ['not_exported_test_status_unknown','NA_pvalue','NA_padj_or_independent_filter','NA_LFC','significant'],default='eligible_nonsignificant')
        verified=table.attrs.get('direction_verified',False)
        rows['direction_verified']=verified
        rows['biological_direction']=rows.source_direction.replace({'positive_LFC':'up','negative_LFC':'down'}) if verified else 'unresolved'
        all_results.append(rows)
        for target in state['rules'].target_group.drop_duplicates():
            subset=rows.loc[rows.target_group.eq(target)]
            sig=subset.loc[subset.significant]
            positive=sig.loc[sig.log2FoldChange.gt(0)];negative=sig.loc[sig.log2FoldChange.lt(0)]
            summaries.append(dict(target_group=target,contrast=name,candidate_genes=len(subset),genes_with_annotation=len(subset),
                                  primary_support_genes=int(subset.primary_support.sum()),in_DE_export=int(subset.in_DE_export.sum()),
                                  tested=int(subset.tested.sum()),eligible=int(subset.eligible.sum()),significant=len(sig),
                                  positive_LFC=len(positive),negative_LFC=len(negative),direction_verified=verified,
                                  biological_up=len(positive) if verified else np.nan,biological_down=len(negative) if verified else np.nan,
                                  not_exported=int((~subset.in_DE_export).sum()),NA_pvalue=int((subset.in_DE_export & subset.pvalue.isna()).sum()),
                                  NA_padj=int((subset.in_DE_export & subset.padj.isna()).sum()),NA_LFC=int((subset.in_DE_export & subset.log2FoldChange.isna()).sum()),
                                  significant_genes=';'.join(sorted(sig['query'])),positive_LFC_genes=';'.join(sorted(positive['query'])),negative_LFC_genes=';'.join(sorted(negative['query'])),
                                  naming_status=naming_status(target)))
            directory=OUT/'contrast_gene_lists'/slug(target);directory.mkdir(parents=True,exist_ok=True)
            # up/down are explicit aliases for signs in the *supplied* DE table,
            # not a claim about a treatment/reference whose polarity is unknown.
            roles={'significant':sig,'up':positive,'down':negative,'NA_or_not_exported':subset.loc[~subset.eligible]}
            for role,genes in roles.items():
                path=directory/f'{name}__{role}.txt'
                path.write_text(''.join(gene+'\n' for gene in sorted(genes['query'])),encoding='utf-8')
                list_manifest.append(dict(path=path.relative_to(OUT).as_posix(),target_group=target,contrast=name,role=role,genes=len(genes),
                                          direction_basis='up=positive supplied LFC; down=negative supplied LFC; treatment/reference polarity unresolved' if not verified else 'verified numerator/denominator LFC',direction_verified=verified))
    combined=pd.concat(all_results,ignore_index=True) if all_results else pd.DataFrame()
    summary=pd.DataFrame(summaries)
    save_table(combined,OUT/'candidate_DE_results.tsv.gz')
    save_table(summary,OUT/'target_contrast_summary.tsv')
    save_table(pd.DataFrame(list_manifest),OUT/'gene_list_manifest.tsv')
    state.update(de_results=combined,de_summary=summary)
    return summary


def plot_candidates(state):
    cfg,summary=state['cfg'],state['summary']
    target_order=state['rules'].target_group.drop_duplicates().tolist()
    counts=summary.set_index('target_group').reindex(target_order)
    fig,ax=plt.subplots(figsize=(9,12))
    sns.heatmap(counts[['GO_genes','KEGG_pathway_genes','keyword_primary_genes','keyword_OG_genes']],annot=True,fmt='g',cmap='Blues',ax=ax)
    ax.set(title='Candidate support per target and evidence type (overlapping gene sets)',xlabel='Evidence type; gene counts are not additive',ylabel='Target')
    savefig(fig,'target_evidence_counts')
    if not state['de_summary'].empty:
        fig,axes=plt.subplots(1,2,figsize=(17,13),layout='constrained')
        for column,title,ax in [('tested','Candidate genes with finite raw pvalue',axes[0]),('significant','Candidate genes passing padj and LFC thresholds',axes[1])]:
            matrix=state['de_summary'].pivot(index='target_group',columns='contrast',values=column).reindex(index=target_order,columns=list(state['tables']))
            sns.heatmap(matrix,annot=True,fmt='g',cmap='YlGnBu',ax=ax)
            ax.set(title=title,ylabel='Target (overlapping membership)',xlabel='Existing DE source label')
        savefig(fig,'target_contrast_counts')
    overlap=pd.DataFrame(0,index=target_order,columns=target_order)
    for a in target_order:
        for b in target_order:
            overlap.loc[a,b]=len(state['gene_sets'][a]&state['gene_sets'][b])
    save_table(overlap,OUT/'target_overlap_matrix.tsv',index=True)
    fig,ax=plt.subplots(figsize=(13,11))
    sns.heatmap(np.log1p(overlap),cmap='mako',ax=ax,cbar_kws={'label':'log(1 + shared genes)'})
    ax.set_title('Overlapping candidate targets; raw shared counts exported separately')
    savefig(fig,'target_overlap')
    genes=sorted(set(state['membership']['query']))
    lfc=pd.DataFrame({name:table.log2FoldChange.reindex(genes) for name,table in state['tables'].items()},index=pd.Index(genes,name='query'))
    sig=pd.DataFrame({name:(table.eligible & table.padj.lt(cfg['padj_threshold']) & table.log2FoldChange.abs().ge(cfg['lfc_threshold'])).reindex(genes,fill_value=False) for name,table in state['tables'].items()},index=lfc.index)
    save_table(lfc,OUT/'candidate_LFC_matrix.tsv',index=True)
    save_table(sig,OUT/'candidate_significance_matrix.tsv',index=True)
    if not lfc.empty and lfc.shape[1]>0:
        annotation=state['catalog'].set_index('query')
        targets=state['membership'].groupby('query').target_group.agg(lambda x:sorted(set(x)))
        sizes={target:len(value) for target,value in state['gene_sets'].items()}
        info=pd.DataFrame({'query':genes})
        info['target_groups']=info['query'].map(targets.map(';'.join))
        info['display_group']=info['query'].map(targets.map(lambda values:min(values,key=lambda t:(sizes[t],target_order.index(t)))))
        info['Preferred_name']=info['query'].map(annotation.Preferred_name)
        info['max_abs_LFC']=info['query'].map(lfc.abs().max(axis=1))
        info['significant_sources']=info['query'].map(sig.sum(axis=1))
        rank=info.sort_values(['max_abs_LFC','query'],ascending=[False,True],na_position='last',kind='stable')
        limit=cfg['lfc_heatmap_max_rows'];quota=max(1,limit//len(target_order))
        selected=[]
        for target in target_order:
            candidates=rank.loc[rank['query'].isin(state['gene_sets'][target]),'query'].head(quota)
            selected.extend(candidates.tolist())
        selected=list(dict.fromkeys(selected))[:limit]
        selected += [gene for gene in rank['query'] if gene not in set(selected)][:max(0,limit-len(selected))]
        info['displayed']=info['query'].isin(selected)
        info['group_order']=info.display_group.map({target:i for i,target in enumerate(target_order)})
        shown=info.loc[info.displayed].sort_values(['group_order','max_abs_LFC','query'],ascending=[True,False,True],kind='stable').reset_index(drop=True)
        save_table(info.drop(columns='group_order'),OUT/'LFC_heatmap_gene_annotations.tsv')
        save_table(shown,OUT/'displayed_LFC_heatmap_rows.tsv')
        plot=lfc.loc[shown['query']]
        fig=plt.figure(figsize=(11,max(9,len(shown)*.065)))
        grid=fig.add_gridspec(1,2,width_ratios=[.15,5],wspace=.07)
        axc=fig.add_subplot(grid[0,0]);ax=fig.add_subplot(grid[0,1])
        palette=dict(zip(target_order,sns.color_palette('husl',len(target_order))))
        axc.imshow(np.array([palette[target] for target in shown.display_group])[:,None,:],aspect='auto',origin='upper')
        axc.set_xticks([]);axc.set_yticks([])
        ax.set_facecolor('#bcbcbc')
        sns.heatmap(plot,mask=plot.isna(),cmap='vlag',center=0,vmin=-cfg['lfc_heatmap_display_limit'],vmax=cfg['lfc_heatmap_display_limit'],
                    yticklabels=False,ax=ax,cbar_kws={'label':'Supplied log2 fold change (display clipped only)'})
        for row,gene in enumerate(shown['query']):
            for column,name in enumerate(plot):
                if sig.loc[gene,name]:
                    ax.text(column+.5,row+.5,'*',ha='center',va='center',fontsize=6,color='black')
        blocks=[]
        for target,block in shown.groupby('display_group',sort=False):
            first,last=block.index.min(),block.index.max()
            blocks.append((target,(first+last)/2))
            ax.axhline(last+1,color='grey',linewidth=.5)
        label_positions=[center for _,center in blocks]
        for i in range(1,len(label_positions)):
            label_positions[i]=max(label_positions[i],label_positions[i-1]+2)
        if label_positions:
            label_positions[-1]=min(label_positions[-1],len(shown)-1)
        for i in range(len(label_positions)-2,-1,-1):
            label_positions[i]=min(label_positions[i],label_positions[i+1]-2)
        for (target,center),position in zip(blocks,label_positions):
            axc.annotate(target,xy=(-.45,center),xytext=(-.8,position),ha='right',va='center',fontsize=7,
                         arrowprops={'arrowstyle':'-','color':'grey','lw':.4},annotation_clip=False)
        ax.set(ylabel='',xlabel='DE source label; biological sign orientation unresolved where unverified',
               title=f'{len(shown)} displayed / {len(lfc)} unique candidate genes; * significant, grey missing\nRow grouping is for display; all overlapping memberships are exported')
        ax.tick_params(axis='x',rotation=45)
        savefig(fig,'candidate_LFC_heatmap')
        state.update(lfc=lfc,significance=sig,heatmap_info=info)
    return list(FIG.glob('*.png'))


def file_expression(state):
    counts,files=state['counts'],state['files']
    if counts is None or files is None:
        state['expression_status']['file_expression']='blocked'
        return state['expression_status']
    complete=counts.loc[~counts.isna().any(axis=1)]
    totals=complete.sum(axis=0)
    if totals.le(0).any():
        state['expression_status']['file_expression']='blocked: nonpositive library total'
        return state['expression_status']
    genes=sorted(set(state['membership']['query'])&set(complete.index))
    keys=files.sort_values(['EXP','ID'],kind='stable').file_key.tolist()
    transformed=np.log2(complete.loc[genes,keys].div(totals.loc[keys],axis=1)*1e6+1)
    transformed.index.name='query'
    std=transformed.std(axis=1,ddof=0)
    scaled=transformed.sub(transformed.mean(axis=1),axis=0).div(std.replace(0,np.nan),axis=0).fillna(0)
    save_table(transformed,OUT/'candidate_file_log2CPM.tsv.gz',index=True)
    save_table(scaled,OUT/'candidate_file_zscores.tsv.gz',index=True)
    variable=transformed.var(axis=1).sort_values(ascending=False,kind='stable')
    selected=variable.loc[variable.gt(0)].head(state['cfg']['file_expression_max_genes']).index
    save_table(pd.DataFrame({'query':selected}),OUT/'displayed_file_expression_genes.tsv')
    methods=dict(transform='log2(expected_count / total_expected_count * 1e6 + 1); z-score each gene across files, population SD',
                 count_type='fractional RSEM gene-level expected_count; no rounding; TPM/FPKM unused',
                 interpretation='file-level exploratory abundance; not a fitted VST/rlog or a count-model analysis',
                 filtering='complete shared gene rows; library denominator includes all complete genes, not just candidates',
                 constant_candidate_genes=int(std.eq(0).sum()),full_candidate_genes=len(transformed),displayed_genes=len(selected),
                 labels='EXP directory membership and file IDs only; treatment/batch sample labels are blocked by missing internal sample identity')
    write_json(OUT/'file_expression_methods.json',methods)
    if len(selected):
        fig,(top,ax)=plt.subplots(2,1,figsize=(15,10),gridspec_kw={'height_ratios':[.25,8]},sharex=True)
        experiments=files.set_index('file_key').loc[keys,'EXP']
        palette=dict(zip(sorted(experiments.unique()),sns.color_palette('Set2',experiments.nunique())))
        top.imshow(np.array([palette[value] for value in experiments])[None,:,:],aspect='auto',extent=(0,len(keys),0,1))
        top.set_yticks([]);top.set_title('File-level expression only: EXP directory colors; sample treatments unverified')
        sns.heatmap(scaled.loc[selected,keys],cmap='vlag',center=0,vmin=-2.5,vmax=2.5,yticklabels=False,xticklabels=[key.split('/')[-1].removesuffix('.genes.results') for key in keys],ax=ax,cbar_kws={'label':'Within-gene z-score of exploratory log2(CPM+1)'})
        ax.tick_params(axis='x',rotation=90,labelsize=6)
        ax.set(ylabel=f'{len(selected)} displayed genes; full candidate matrix exported',xlabel='RSEM file identity (not confirmed biological sample)')
        savefig(fig,'candidate_file_expression')
    state['expression_status']['file_expression']='complete'
    write_json(OUT/'expression_mapping_status.json',state['expression_status'])
    return methods


def verify_outputs(state=None,check_notebook=False):
    cfg=json.loads((OUT/'configuration_used.json').read_text())
    membership=pd.read_csv(OUT/'gene_target_membership.tsv',sep='\t',keep_default_na=False)
    evidence=pd.read_csv(OUT/'gene_match_evidence.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
    summary=pd.read_csv(OUT/'target_gene_counts.tsv',sep='\t')
    de=pd.read_csv(OUT/'candidate_DE_results.tsv.gz',sep='\t')
    lists=pd.read_csv(OUT/'gene_list_manifest.tsv',sep='\t',keep_default_na=False)
    checks=[]
    def check(name,value):
        checks.append(dict(check=name,passed=bool(value)))
    expected=['configuration_used.json','preserved_files_sha256.tsv','input_headers_and_paths.json','contrasts_used.json','input_issues.json',
              'annotation_parser_audit.json','duplicated_annotation_queries.tsv','count_file_audit.tsv','sample_file_mapping.tsv',
              'control_metadata.tsv','metadata_design_counts.tsv','expression_mapping_status.json','candidate_rules.tsv',
              'synonym_expansions_and_patterns.tsv','gene_match_evidence.tsv.gz','gene_target_membership.tsv','target_gene_counts.tsv',
              'target_evidence_type_counts.tsv','matched_gene_catalog.tsv','rule_match_counts.tsv','identifier_mapping_diagnostics.tsv',
              'unmatched_identifiers.tsv','target_overlap_pairs.tsv','target_overlap_matrix.tsv','candidate_DE_results.tsv.gz',
              'target_contrast_summary.tsv','gene_list_manifest.tsv','candidate_LFC_matrix.tsv','candidate_significance_matrix.tsv',
              'LFC_heatmap_gene_annotations.tsv','displayed_LFC_heatmap_rows.tsv','candidate_file_log2CPM.tsv.gz','candidate_file_zscores.tsv.gz',
              'displayed_file_expression_genes.tsv','file_expression_methods.json','synthetic_smoke_tests.json','synthetic_smoke_tests.txt']
    for name in expected:
        check(f'file:{name}',(OUT/name).is_file() and (OUT/name).stat().st_size>0)
    for name in ['target_evidence_counts','target_contrast_counts','target_overlap','candidate_LFC_heatmap','candidate_file_expression']:
        for extension in ['png','pdf']:
            check(f'figure:{name}.{extension}',(FIG/f'{name}.{extension}').is_file() and (FIG/f'{name}.{extension}').stat().st_size>0)
    check('unique_query_target_memberships',not membership.duplicated(['query','target_group']).any())
    check('evidence_routes_preserved',set(zip(membership['query'],membership.target_group))==set(zip(evidence['query'],evidence.target_group)))
    check('no_KO_assignment_used_as_pathway_evidence',not evidence.matched_field.eq('KEGG_ko').any())
    for kind,field in [('GO','GOs'),('KEGG_reference_pathway','KEGG_Pathway')]:
        subset=evidence.loc[evidence.annotation_type.eq(kind)]
        check(f'exact_membership:{kind}',subset.matched_field.eq(field).all() and all(r.effective_alias in r.matched_value.split(',') for r in subset.itertuples()))
    keyword=evidence.loc[evidence.annotation_type.eq('gene_family_keyword')]
    check('keywords_only_allowed_fields',keyword.matched_field.isin(['Preferred_name','Description','eggNOG_OGs']).all())
    check('OG_evidence_marked_supplementary',keyword.loc[keyword.matched_field.eq('eggNOG_OGs'),'evidence_tier'].eq('supplementary_OG_only_route').all())
    check('Eco_remains_unresolved',summary.loc[summary.target_group.eq('Eco'),'naming_status'].eq('unresolved_literal_query').all())
    check('HyAlx_requires_confirmation',summary.loc[summary.target_group.eq('HyAlx'),'naming_status'].eq('project_specific_confirmation_required').all())
    actual_headers=json.loads((OUT/'input_headers_and_paths.json').read_text())
    check('actual_candidate_header_preserved',actual_headers['candidate']==list(pd.read_csv(ROOT/cfg['candidate_file'],nrows=0)))
    for row in summary.itertuples():
        check(f'target_count:{row.target_group}',row.genes==membership.loc[membership.target_group.eq(row.target_group),'query'].nunique())
    source_items=json.loads((OUT/'contrasts_used.json').read_text())
    for item in source_items:
        name=item['name']
        source=read_contrast(item)
        check(f'actual_DE_header_preserved:{name}',actual_headers['DE'][name]==list(pd.read_csv(ROOT/item['path'],sep=item.get('delimiter'),engine='python',nrows=0)))
        subset=de.loc[de.contrast.eq(name)].drop_duplicates('query').set_index('query')
        check(f'candidate_DE_unique_per_target:{name}',not de.loc[de.contrast.eq(name)].duplicated(['query','target_group']).any())
        for column in ['log2FoldChange','pvalue','padj']:
            expected_values=source[column].reindex(subset.index).to_numpy(dtype=float)
            check(f'DE_statistics_preserved:{name}/{column}',np.allclose(subset[column].to_numpy(dtype=float),expected_values,rtol=1e-8,atol=1e-12,equal_nan=True))
        check(f'NA_not_nonsignificant:{name}',subset.loc[~subset.eligible,'DE_status'].ne('eligible_nonsignificant').all())
    for row in lists.itertuples():
        path=OUT/row.path
        check(f'gene_list_exists:{row.path}',path.is_file())
        if not path.is_file():
            continue
        genes=path.read_text(encoding='utf-8').splitlines()
        if not row.contrast:
            wanted=set(membership.loc[membership.target_group.eq(row.target_group),'query'])
        else:
            candidates=de.loc[de.target_group.eq(row.target_group) & de.contrast.eq(row.contrast)]
            mask={'significant':candidates.significant,'up':candidates.significant & candidates.log2FoldChange.gt(0),
                  'down':candidates.significant & candidates.log2FoldChange.lt(0),'NA_or_not_exported':~candidates.eligible}[row.role]
            wanted=set(candidates.loc[mask,'query'])
        check(f'gene_list_content:{row.path}',set(genes)==wanted and len(genes)==row.genes and genes==sorted(set(genes)))
    baseline=pd.read_csv(OUT/'preserved_files_sha256.tsv',sep='\t')
    for row in baseline.itertuples():
        check(f'preserved:{row.path}',(ROOT/row.path).is_file() and sha256(ROOT/row.path)==row.sha256)
    check('protected_inventory_unchanged',{p.relative_to(ROOT).as_posix() for p in protected_paths()}==set(baseline.path))
    smoke=json.loads((OUT/'synthetic_smoke_tests.json').read_text())
    check('synthetic_matching_smoke_tests',smoke['passed'])
    if check_notebook:
        document=json.loads(NOTEBOOK.read_text(encoding='utf-8'))
        cells=[cell for cell in document['cells'] if cell['cell_type']=='code']
        counters=[cell.get('execution_count') for cell in cells]
        errors=[output for cell in cells for output in cell.get('outputs',[]) if output.get('output_type')=='error']
        check('clean_notebook_execution',bool(cells) and counters==list(range(1,len(cells)+1)) and not errors)
        write_json(OUT/'notebook_execution.json',dict(path=NOTEBOOK.relative_to(ROOT).as_posix(),execution_counts=counters,errors=len(errors),sha256=sha256(NOTEBOOK)))
    table=pd.DataFrame(checks)
    save_table(table,OUT/'output_validation.tsv')
    result=dict(checks=len(checks),passed=bool(table.passed.all()),failed=table.loc[~table.passed,'check'].tolist(),notebook_checked=check_notebook)
    write_json(OUT/'output_validation.json',result)
    if not result['passed']:
        raise AssertionError(result)
    return result


def write_report():
    cfg=json.loads((OUT/'configuration_used.json').read_text())
    summary=pd.read_csv(OUT/'target_gene_counts.tsv',sep='\t')
    membership=pd.read_csv(OUT/'gene_target_membership.tsv',sep='\t',keep_default_na=False)
    results=pd.read_csv(OUT/'candidate_DE_results.tsv.gz',sep='\t')
    validation=json.loads((OUT/'output_validation.json').read_text())
    audit=json.loads((OUT/'annotation_parser_audit.json').read_text())
    expression=json.loads((OUT/'expression_mapping_status.json').read_text())
    source_items=json.loads((OUT/'contrasts_used.json').read_text())
    rules=pd.read_csv(OUT/'candidate_rules.tsv',sep='\t')
    issues=json.loads((OUT/'input_issues.json').read_text())
    packages={name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy','matplotlib','seaborn']}
    write_json(OUT/'environment_versions.json',dict(python=sys.version,executable=sys.executable,packages=packages))
    lines=['# Development pathways and gene-family candidate analysis','',
           f"{len(rules)} supplied rules / {summary.target_group.nunique()} target groups yield **{membership['query'].nunique()} distinct candidate genes** and {len(membership)} overlapping gene/target memberships. Candidate support does not establish pathway activity, sequence orthology, or signaling activation.",'',
           '## Inputs and matching','',
           f"Candidate source: `{cfg['candidate_file']}`; guidance: `{cfg['guidance_file']}`. The exact user-specified files exist; no candidate-file substitution occurred. Annotation source: `{cfg['annotation_file']}`. The task's `results/deseq2_rsem/tables/` directory is absent. The verified existing DE inputs are in `data/deseq2_tables/tables/`, as recorded by the pre-existing contrast manifest; their headers and full contents were validated before analysis.",'',
           f"The eggNOG header contains query, Description, Preferred_name, GOs, KEGG_ko, KEGG_Pathway and eggNOG_OGs. There are {audit['rows']} records / {audit['unique_queries']} unique queries / {audit['duplicate_queries']} duplicate queries. Comments, #query, tab-delimited fields, comma-delimited terms and missing '-' are parsed explicitly. Duplicate records retain record numbers in the evidence table; memberships union exact query/target pairs. Full headers are in `input_headers_and_paths.json`.",'',
           'GO rules match exact comma-token IDs in GOs. KEGG reference rules match exact koNNNNN/mapNNNNN representations of the same five-digit pathway code in KEGG_Pathway. Both representations are retained as evidence but each gene/target is counted once. KEGG_ko is reported as context only: individual KO assignments such as ko:K07188 are never searched as pathway maps. No pathway-to-K-number expansion or web-derived synonym expansion was performed. Supplied source URLs are retained as provenance, not treated as retrieved ontology definitions.','',
           'Keyword rules use only pipe-separated aliases supplied by the candidate table, case-insensitively in Preferred_name/Description. eggNOG_OGs is supplementary only. Token/phrase boundaries suppress accidental substring hits (e.g. Eco inside ecology). Listed FOX/HOX stems permit a family-letter/numeric symbol suffix; SOX/ALX/ZIC/TDRD/NANOS/PIWI/GSC stems permit numeric/letter suffixes to recover explicitly named family members. These matching patterns and effective aliases are exported in `synonym_expansions_and_patterns.tsv`. Case-equivalent alias duplicates collapse; no new biological synonyms are invented. Generic homeobox/T-box/forkhead/Argonaute matches remain broad family candidates, not one-to-one ortholog assignments.','',
           '**Eco remains unresolved**: only literal Eco-token evidence is retained, with no assumed relationship to EcoKD1. **HyAlx requires project-specific identity/orthology confirmation**: only supplied HyAlx/Alx/ALX aliases are used; generic homeobox alone is not added to HyAlx. All candidate memberships are marked for manual review.','',
           'Zero hits under a supplied rule do not establish biological absence. Annotations can be incomplete, literal aliases can be missing, and project-specific naming or sequence-level orthology work may be necessary.','',
           '## Target counts by evidence','',
           'Counts are unique genes within each target/evidence category and are not additive. GO, KEGG and keyword evidence can support the same gene; target groups remain separate. The supplementary OG column includes any OG-route support; OG-only genes have no primary matching field.','',
           '| Target | All genes | GO | KEGG pathway | Keyword primary | Supplementary OG | OG only |','|---|---:|---:|---:|---:|---:|---:|']
    for row in summary.itertuples():
        lines.append(f'| {row.target_group} | {row.genes} | {row.GO_genes} | {row.KEGG_pathway_genes} | {row.keyword_primary_genes} | {row.keyword_OG_genes} | {row.supplementary_only_genes} |')
    lines += ['', 'BMP, TGF-beta, TGF-beta superfamily, SMAD and TGF-beta/BMP remain distinct. FGF, MAPK, JNK and MAPK/FGF/JNK also remain distinct. `target_overlap_pairs.tsv` lists every shared gene and Jaccard score; `target_overlap_matrix.tsv` contains shared-gene counts. Unique overall totals deduplicate queries across targets.','',
              '## Existing DESeq2 statistics and key results','',
              f"Thresholds: padj < {cfg['padj_threshold']}, |log2FoldChange| >= {cfg['lfc_threshold']}, with finite raw pvalue, padj and LFC. Statistics are reused verbatim from the original CSVs; no R/PyDESeq2 fitting, rounding, or count-model inference occurs. Genes missing from an exported table have unknown test status. NA pvalues, independent-filtered/undefined padj and missing LFCs remain separate from eligible nonsignificant genes.",'',
              'The existing source labels are EcoKD1, 3OHC12, Conventionalized, Cvbct and Wild. The files contain statistics but no embedded numerator/denominator or model design. Their direction_verified flags remain false. Positive/negative supplied-LFC counts are valid numerical summaries; **biological treatment-direction interpretation is blocked**. The requested `__up.txt`/`__down.txt` files are explicit aliases for positive/negative significant LFCs in the supplied source table only, not claims about a confirmed treatment/reference. Every list records this basis in `gene_list_manifest.tsv`; biological_up/down summary fields are undefined where polarity is unverified.','',
              '| Source | Distinct candidate genes tested | Eligible | Significant | Positive LFC | Negative LFC |','|---|---:|---:|---:|---:|---:|']
    drivers=[]
    for item in source_items:
        name=item['name'];subset=results.loc[results.contrast.eq(name)].drop_duplicates('query')
        sig=subset.loc[subset.significant]
        lines.append(f'| {name} | {int(subset.tested.sum())} | {int(subset.eligible.sum())} | {int(subset.significant.sum())} | {int((sig.log2FoldChange>0).sum())} | {int((sig.log2FoldChange<0).sum())} |')
        top=sig.sort_values(['padj','query'],kind='stable').head(3)
        for row in top.itertuples():
            drivers.append(f"- {name}: `{row.query}` ({row.Preferred_name if isinstance(row.Preferred_name,str) and row.Preferred_name else 'unnamed'}), LFC={row.log2FoldChange:.3g}, padj={row.padj:.3g}.")
    lines += ['', 'Selected statistical drivers (overlapping target assignments are retained in the full evidence and DE tables):','',*drivers,'',
              'Per-target/per-contrast counts and exact driving-gene lists are in `target_contrast_summary.tsv`; `candidate_DE_results.tsv.gz` retains every candidate membership, original statistics, source-sign direction and explicit NA/export status. Readable symbols/descriptions and all original annotation context are in those tables and the per-target TSV lists.','',
              '## Expression and plots','',
              f"RSEM inputs contain gene_id, transcript_id(s), length, effective_length, expected_count, TPM and FPKM. All count genes/sample-file keys are validated for uniqueness and numeric validity. The current 61 files share 28,917 gene IDs and contain 194,391 fractional expected-count values. No fractions are rounded; TPM and FPKM are unused. Identity joins use exact HVAEP<scaffold>.G identifiers; no row-order joins, transcript substitution or suffix stripping occur. Unmatched IDs and parser duplicates are explicitly exported.",'',
              f"**Treatment-labelled sample expression is blocked**: {expression['reason']} There are {expression['filename_candidate_matches']} candidate metadata/file matches, {expression['files_without_metadata']} files without metadata and {expression['metadata_without_files']} metadata records without a file. Metadata EXP directories and treatment/control category membership are verified. Filename candidates remain unverified biological identities. No samples are automatically removed.",'',
              'The independent file-expression view uses log2(expected_count / total_expected_count * 1e6 + 1), with the denominator from all complete shared genes, then scales each gene across files by its population SD. Constant genes receive zero z-scores in exports and are not chosen for the variance-ranked display. This is an exploratory compositional transform, not a fitted VST/rlog; residual mean-variance dependence remains. Files are labelled by their IDs and EXP directory membership only. Full candidate transformed/z-score matrices and displayed rows are exported.','',
              'The LFC heatmap exports the full matched-gene matrix and significance mask. Its displayed subset is deterministic: target-balanced high-absolute-LFC candidates followed by remaining highest-absolute-LFC genes. Unique gene rows are grouped by their smallest supporting target for display only; all overlapping targets are retained in the row annotation table. Stars mark cells passing both thresholds; grey is missing. Color clipping affects display only. Target/evidence counts, tested/significant counts and target-overlap plots use full sets.','',
              'Optional candidate clustering and enrichment were not performed: neither is necessary for the requested evidence extraction/DE summaries, and avoiding extra pathway-activity inference is appropriate for these transferred predictions.','',
              '## Validation, preservation and files','',
              'Synthetic tests cover actual comment/header/missing-data parsing, duplicated queries, exact GO tokens, equivalent KEGG map forms, refusal to match K04310 as pathway ko04310, family symbol suffixes, keyword substring false positives, OG-only evidence, Eco/HyAlx flags, row-order invariance and shuffled DE/NA joins.','',
              f"Validation: **{validation['checks']} checks; passed={validation['passed']}**; notebook_checked={validation['notebook_checked']}. `output_validation.tsv` includes all expected outputs, every gene-list content/count, source-statistic comparisons, namespace checks, and SHA-256 preservation checks for all raw inputs and pre-existing scripts/notebooks/results. Existing analysis outputs were not regenerated or modified. No dependency change was required.",'',
              'Notebook: `scripts/05_development_pathway_analysis.ipynb`. New reusable modules/config are under scripts/. All candidate tables, evidence, diagnostics, gene lists, manifests and this report are under `results/processedData/development_pathways/`; PNG/PDF figures are under `results/figures/development_pathways/`. `output_inventory.tsv` lists each generated artifact and byte size; empty lists for zero-hit targets are intentional and validated.','',
              f'Additional input issues: {issues if issues else "none beyond the explicitly documented sample identity and contrast-direction limitations"}.','']
    (OUT/'analysis_report.md').write_text('\n'.join(lines),encoding='utf-8')
    records=[]
    for directory in [OUT,FIG]:
        for path in sorted(directory.rglob('*')):
            if path.is_file() and path.name!='output_inventory.tsv':
                records.append(dict(path=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size))
    for name in sorted(NEW_SCRIPTS):
        path=ROOT/'scripts'/name
        if path.is_file():
            records.append(dict(path=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size))
    save_table(pd.DataFrame(records),OUT/'output_inventory.tsv')
    return OUT/'analysis_report.md'


def run_analysis(overrides=None):
    from development_pathway_smoke_tests import run_tests
    state=start_analysis(overrides)
    extract_candidates(state)
    analyze_de(state)
    plot_candidates(state)
    file_expression(state)
    run_tests()
    verify_outputs(state)
    print('Report:',write_report())
    return state


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-only',action='store_true',help='Verify saved artifacts and clean notebook execution without recomputing analyses')
    args=parser.parse_args()
    if args.verify_only:
        print(verify_outputs(check_notebook=True))
        print('Report:',write_report())
    else:
        run_analysis()
