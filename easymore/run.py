import configparser
import argparse
from easymore_parallel import esmr

parser = argparse.ArgumentParser(description='Run easymore')
parser.add_argument('-c', metavar='c', type=str, required=True,
                    help='configuration file')
args= parser.parse_args()

config= configparser.ConfigParser()
config.read(args.c)
inputs= config['input']

esmr(
    case_name=inputs['case_name'],
    target_shp= inputs['target_shp'],
    target_shp_ID=inputs['target_shp_ID'],
    target_shp_lat=inputs['target_shp_lat'],
    target_shp_lon=inputs['target_shp_lon'],
    source_nc=inputs['source_nc'],
    var_names= inputs['var_names'],
    var_lon= inputs['var_lon'],
    var_lat= inputs['var_lat'],
    var_time= inputs['var_time'],
    var_ID= inputs['var_ID'],
    var_names_remapped= inputs['var_names_remapped'],
    source_shp= inputs['source_shp'],
    source_shp_lat= inputs['source_shp_lat'],
    source_shp_lon= inputs['source_shp_lon'],
    source_shp_ID= inputs['source_shp_ID'],
    remapped_var_id= inputs['remapped_var_id'],
    remapped_var_lat= inputs['remapped_var_lat'],
    remapped_var_lon= inputs['remapped_var_lon'],
    remapped_dim_id= inputs['remapped_dim_id'],
    overwrite_existing_remap= inputs['overwrite_existing_remap'],
    temp_dir= inputs['temp_dir'],
    output_dir= inputs['output_dir'],
    format_list=inputs['format_list'],
    fill_value_list= inputs['fill_value_list'],
    remap_csv= inputs['remap_csv'],
    author_name= inputs['author_name'],
    license= inputs['license'],
    save_csv= inputs['save_csv'],
    save_shp= inputs['save_shp'],
    ncpus= inputs['ncpus'],
    sort_ID= inputs['sort_ID'],
    complevel= inputs['complevel']
)
