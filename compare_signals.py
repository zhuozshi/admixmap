import subprocess
import pandas as pd
import numpy as np
import argparse
import sys

def pull_variants_ld(variant_list: list,anc: str, ld_link_token: str, output_path: str, genome_build='grch38_high_coverage'):
    if output_path[-1]!='/':
        output_path+='/'
    for j in variant_list:
        command = (f"curl -k -X GET 'https://ldlink.nih.gov/LDlinkRest/ldproxy?var={j}&pop={anc}&r2_d=r2&window=1000000&genome_build={genome_build}&token={ld_link_token}' > {output_path}ldlink_{j}.{anc}.txt")
        subprocess.run(command, shell=True, check=True)

def read_ld_variants(variant_list: list, input_path: str, ancs='EUR',r2_cutoff=0.8):
    variants_notin_ref = []
    inld = ''
    if ancs == 'EACH':
        ancs = ['AFR','AMR','SAS','EAS','EUR','ALL']
    else:
        ancs = [ancs]
    if input_path[-1]!='/':
        input_path+='/'
    if r2_cutoff>1 or r2_cutoff<0:
        print('ERROR: R2 threshold must be between 0 and 1')
        raise()
    for j in variant_list:
        for anc in ancs:
            ld_proxy_out = pd.read_table(f'{input_path}ldlink_{j}.{anc}.txt')
            if 'RS_Number' in ld_proxy_out.columns:
                ld_proxy_out = ld_proxy_out[ld_proxy_out['R2']>r2_cutoff]
                if len(ld_proxy_out)>0:
                    ld_proxy_out['lead_variant'] = j
                    ld_proxy_out['ancestry'] = anc
                    if type(inld) == str:
                        inld = ld_proxy_out.copy()
                    else:
                        inld = pd.concat([inld,ld_proxy_out])
                else:
                    variants_notin_ref.append(j)
            else:
                variants_notin_ref.append(j)
    variants_notin_ref = set(variants_notin_ref)
    return(inld,variants_notin_ref)

def read_reference_assocs(ref_assoc_path: str, chrom_colname, pos_colname, sep=','):
    ref_assoc = pd.read_table(ref_assoc_path,sep=sep)
    column_header_inputs = [chrom_colname,pos_colname]
    for i in column_header_inputs:
        if i not in ref_assoc.columns:
            print(f'ERROR: {i} could not be found within {ref_assoc_path}')
            raise()

    if type(ref_assoc[chrom_colname].to_list()[0]) != str:
        ref_assoc[chrom_colname] = 'chr'+ref_assoc[chrom_colname].astype(str)
    if 'chr' not in ref_assoc[chrom_colname].to_list()[0]:
        ref_assoc[chrom_colname] = 'chr'+ref_assoc[chrom_colname].astype(str)
    ref_assoc['id'] = ref_assoc[chrom_colname]+':'+ref_assoc[pos_colname].astype(str)
    return(ref_assoc)

def read_new_assocs(new_assoc_path: str, chrom_colname, start_colname, id_colname, end_colname=None, sep=','):
    # two input types are possible. 1) a dataframe of variants 2) a dataframe of regions (i.e. peaks, loci...)
    new_assoc = pd.read_table(new_assoc_path,sep=sep)
    column_header_inputs = [chrom_colname,start_colname,id_colname]
    if end_colname != None:
        column_header_inputs.append(end_colname)
    for i in column_header_inputs:
        if i not in new_assoc.columns:
            print(f'ERROR: {i} could not be found within {new_assoc_path}')
            raise()

    if end_colname == None:
        end_colname='end'
        new_assoc[end_colname] = new_assoc[start_colname]

    if type(new_assoc[chrom_colname].to_list()[0]) != str:
        new_assoc[chrom_colname] = 'chr'+new_assoc[chrom_colname].astype(str)
    if 'chr' not in new_assoc[chrom_colname].to_list()[0]:
        new_assoc[chrom_colname] = 'chr'+new_assoc[chrom_colname]

    new_assoc = new_assoc[[chrom_colname,start_colname,end_colname,id_colname]]
    new_assoc.columns = ['chr','start','end','id']
    return(new_assoc)

def compare_signals(new_assoc,variants_notin_ref,inld):
    signal_status = {'id':[],'known-or-novel':[],'known_var':[],'r2':[],'MAF':[]}
    if type(inld) == str:
        return(new_assoc[['id']])
    inld_id = inld['Coord'].str.split(':',expand=True)
    inld['chr'] = inld_id[0]
    inld['pos'] = inld_id[1].astype(int)
    inld = inld[['chr','pos','Coord','lead_variant','ancestry','R2','MAF']]
    vardf = {'chr':[],'pos':[],'Coord':[],'lead_variant':[],'ancestry':[],'R2':[],'MAF':[]}
    for i in variants_notin_ref:
        if i not in inld['lead_variant']:
            vardf['chr'].append(i.split(':')[0])
            vardf['pos'].append(int(i.split(':')[1]))
            vardf['Coord'].append(i)
            vardf['lead_variant'].append(i)
            vardf['ancestry'].append('')
            vardf['R2'].append(1)
            vardf['MAF'].append(np.nan)
    if len(vardf['chr'])>0:
        vardf = pd.DataFrame(vardf)
        inld = pd.concat([inld,vardf])

    for index,row in new_assoc.iterrows():
        t_inld = inld[(inld['chr']==row['chr']) & (inld['pos']>=row['start']) & (inld['pos']<=row['end'])]
        signal_status['id'].append(row['id'])
        if len(t_inld)>0:
            kn = 'known'
            known_var = t_inld['lead_variant'].to_list()[0]
            r2_val = ''
            MAF = ''
            for anc in t_inld['ancestry'].unique():
                tanc_inld = t_inld[t_inld['ancestry']==anc]
                tanc_inld = tanc_inld.sort_values('R2',ascending=False)
                r2_val += f"{anc}: {np.round(tanc_inld['R2'].to_list()[0],2)},"
                if not np.isnan(tanc_inld['MAF'].to_list()[0]):
                    MAF += f"{anc}: {np.round(tanc_inld['MAF'].to_list()[0],3)},"
            r2_val = r2_val[:-1]
            MAF = MAF[:-1]
        else:
            kn = 'novel'
            r2_val = np.nan
            MAF = np.nan
            known_var = np.nan
        signal_status['known-or-novel'].append(kn)
        signal_status['r2'].append(r2_val)
        signal_status['MAF'].append(MAF)
        signal_status['known_var'].append(known_var)

    return(pd.DataFrame(signal_status))

def export_results(signal_status, output: str, sep=','):
    signal_status.to_csv(output,index=None, sep=sep)

def create_reference(ref_assoc_path: str, chrom_colname, pos_colname, output_path: str, ld_link_token: str, var_file='variants.txt',anc='EUR', genome_build='grch38_high_coverage', sep=','):
    ref_assoc = read_reference_assocs(ref_assoc_path=ref_assoc_path, sep=',', chrom_colname=chrom_colname, pos_colname=pos_colname)
    if anc=='EACH':
        ancs = ['EUR','SAS','AFR','AMR','EAS','ALL']
    else:
        ancs = [anc]
    for anc in ancs:
        pull_variants_ld(variant_list=ref_assoc['id'].to_list(),anc=anc,output_path=output_path,ld_link_token=ld_link_token,genome_build=genome_build)
    ref_assoc[['id']].to_csv(f'{output_path}{var_file}',index=None)

def run_comparison(new_assoc_path: str, ld_path: str, output: str, var_file='variants.txt', assoc_chrom_colname='CHR', assoc_start_colname='POS', assoc_id_colname='ID', assoc_end_colname=None, ancs='EUR', r2_cutoff=0.8,assoc_sep=','):
    if ld_path[-1] != '/':
        ld_path += '/'
    new_assoc = read_new_assocs(new_assoc_path, sep=assoc_sep, chrom_colname=assoc_chrom_colname, start_colname=assoc_start_colname, id_colname=assoc_id_colname, end_colname=assoc_end_colname)
    ref_assoc = pd.read_table(f'{ld_path}{var_file}')
    inld, variants_notin_ref = read_ld_variants(ref_assoc['id'].to_list(), input_path=ld_path, ancs=ancs,r2_cutoff=r2_cutoff)
    comparison = compare_signals(new_assoc,variants_notin_ref,inld)
    export_results(comparison, output, sep=',')

def main():

    parser = argparse.ArgumentParser()
    #subparser = parser.add_subparsers()

    subparsers = parser.add_subparsers(title='Commands', dest='command')
    # use dispatch pattern to invoke method with same name

    #parser.add_argument('pull_ld')
    ldpull = subparsers.add_parser('pull_ld')
    ldpull.add_argument("-r", "--reference",dest='ref_assoc_path')
    ldpull.add_argument('-rs','--reference_sep',dest='ref_sep',default=',')
    ldpull.add_argument('-rc','--reference_chrom',dest='ref_chrom_colname',default='CHR')
    ldpull.add_argument('-rp','--reference_pos',dest='ref_pos_colname',default='POS')
    ldpull.add_argument('-g','--genome_build',dest='genome_build',default='grch38_high_coverage')
    ldpull.add_argument('-t','--ldlink_token',dest='token')
    ldpull.add_argument('-anc','--ancestries',dest='ancs',default='EUR')
    ldpull.add_argument('-vi','--reference_var_output',dest='ref_var_out',default='variants.txt')
    ldpull.add_argument('-lo','--ld_output',dest='ld_output_path')

    #parser.add_argument('-c','--compare',dest='compare')
    compare = subparsers.add_parser('compare')
    compare.add_argument("-vi", "--reference_variants",dest='ref_assoc_variants',default='variants.txt')
    compare.add_argument('-a', '--association_file',dest='new_assoc_path')
    compare.add_argument('-as','--association_sep',dest='assoc_sep',default=',')
    compare.add_argument('-ac','--association_chrom',dest='assoc_chrom_colname',default='CHR')
    compare.add_argument('-ap','--association_startpos',dest='assoc_start_colname',default='POS')
    compare.add_argument('-ai','--association_id',dest='assoc_id_colname',default='ID')
    compare.add_argument('-ae','--association_endpos',dest='assoc_end_colname')
    compare.add_argument('-li','--ld_path',dest='ld_path')
    compare.add_argument('-r2','--r2',dest='r2_cutoff',default=0.8)
    compare.add_argument('-anc','--ancestries',dest='ancs',default='EUR')
    compare.add_argument("-o", "--output",dest = "output")


    args = parser.parse_args()
    if args.command == 'pull_ld':
        print('Extracting LD information')
        create_reference(var_file=args.ref_var_out, ref_assoc_path=args.ref_assoc_path, sep=args.ref_sep, chrom_colname=args.ref_chrom_colname, pos_colname=args.ref_pos_colname, genome_build=args.genome_build, output_path=args.ld_output_path, ld_link_token=args.token,anc=args.ancs)
        # Implement logic for the 'foo' command here
    elif args.command == 'compare':
        print('Comparing association results with LD reference')
        run_comparison(var_file=args.ref_assoc_variants, new_assoc_path=args.new_assoc_path, assoc_sep=args.assoc_sep, assoc_chrom_colname=args.assoc_chrom_colname, assoc_start_colname=args.assoc_start_colname, assoc_id_colname=args.assoc_id_colname, assoc_end_colname=args.assoc_end_colname, ld_path=args.ld_path, r2_cutoff=args.r2_cutoff,output=args.output,ancs=args.ancs)
        # Implement logic for the 'bar' command here
    else:
        # If an unknown command is provided, show the help message
        parser.print_help()
    #if args.pull_ld is not None:
        #create_reference(ref_assoc_path=args.ref_assoc_path, sep=args.ref_sep, chrom_colname=args.ref_chrom_colname, pos_colname=args.ref_pos_colname, genome_build=args.genome_build, output_path=args.ld_output_path, ld_link_token=args.token,anc=args.ancs)
    #elif args.compare is not None:
        #run_comparison(new_assoc_path=args.new_assoc_path, assoc_sep=args.assoc_sep, assoc_chrom_colname=args.assoc_chrom_colname, assoc_start_colname=args.assoc_start_colname, assoc_id_colname=args.assoc_id_colname, assoc_end_colname=args.assoc_end_colname, ref_assoc_path=args.ref_assoc_path, ref_sep=args.ref_assoc_path, ref_chrom_colname=args.ref_chrom_colname, ref_pos_colname=args.ref_pos_colname, ld_path=args.ld_path, r2_cutoff=args.r2_cutoff,output=args.output,ancs=args.ancs)
if __name__ == "__main__":
    main()
