#!/usr/bin/env python2
# This work is licensed under the terms of the MIT license. See LICENSE for the full text.

import os
import sys
import re
import shutil
import pandas as pd

if __name__ == "__main__":
    tables_output_dir = os.path.join('analysis_output', 'tables')
    if os.path.isdir(tables_output_dir):
        shutil.rmtree(tables_output_dir)

    analysis_dirs = [os.path.join('analysis_output', d) for d in os.listdir('analysis_output') if os.path.isdir(os.path.join('analysis_output', d))]

    stats_files = [os.path.join(d, '%s-stats.txt' % os.path.basename(d)) for d in analysis_dirs]

    data_dicts = {}
    for fpath in stats_files:
        if not os.path.isfile(fpath):
            print fpath
            raise Exception()

        table_name = None
        new_data = []
        with open(fpath, 'r') as f:
            run_name = os.path.basename(fpath)[:-10]
            for line in f:
                m = re.match('(\S+)(?:\s+[-]\s+)(\S+)(?:\s+vs\s+)(\S+)(?:.*?)', line)
                if m:
                    table_name = '%s-%s-vs-%s.csv' % ( m.group(1), m.group(2), m.group(3) )
                else:
                    if ':' in line:
                        assert( table_name )
                        datum_name, datum_value = line.split(':')
                        datum_name = datum_name.strip()
                        datum_value = float( datum_value.split()[0] )
                        new_data.append( (datum_name, datum_value) )
                    elif len(new_data) > 0:
                        if table_name not in data_dicts:
                            data_dicts[table_name] = {datum_name : [datum_value] for datum_name, datum_value in new_data}
                            data_dicts[table_name]['run_names'] = [run_name]
                            data_dicts[table_name]['table_name'] = [table_name[:-4]]
                        else:
                            assert( run_name not in data_dicts[table_name] )
                            data_dicts[table_name]['run_names'].append(run_name)
                            data_dicts[table_name]['table_name'].append(table_name[:-4])
                            for datum_name, datum_value in new_data:
                                data_dicts[table_name][datum_name].append(datum_value)
                        new_data = []

    os.makedirs(tables_output_dir)

    all_ps = []
    for table_name in data_dicts:
        special_columns = ['run_names', 'table_name', 'n', 'Fraction correct', 'MAE', "Pearson's R", "Spearman's R"]
        all_columns = [x for x in data_dicts[table_name].keys() if x not in special_columns]
        special_columns.extend(all_columns)
        p = pd.DataFrame(data_dicts[table_name], columns=special_columns)
        all_ps.append(p)
        individual_p = p.drop('table_name', axis=1)
        individual_p.set_index('run_names', inplace=True)
        individual_p.to_csv(os.path.join(tables_output_dir, table_name))

    # Make all data table    
    all_df = pd.concat(all_ps)
    all_df.to_csv( os.path.join(tables_output_dir, 'all-data.csv') )


