#!/usr/bin/env python3

# This file is part of minMutFinder.
#
# minMutFinder is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# minMutFinder is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with minMutFinder. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2024 Respiratory Viruses Unit, Microbiology Department, Vall d’Hebron Hospital Universitari, Vall d’Hebron Institut de Recerca (VHIR), Vall d’Hebron Barcelona Hospital Campus, Passeig Vall d’Hebron 119-129, 08035 Barcelona, Spain

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

def visualizer(tsv_check, csv, csv_muts, outd_p, outd_csv, proteins, annotate, SAMPLE, syn_muts, lengths_df):
    proteins_list = proteins
    rows = len(proteins_list)

    if tsv_check:
        qc = pd.read_csv(csv, sep=';')
        qc.columns = ["ref", "test", "Count"]

        qc_T = pd.DataFrame(columns=['protein', 'Length', 'N_percentage', 'Coverage', 'Q1_Depth', 'Median_Depth', 'Q3_Depth'])

        for g in proteins_list:
            df = pd.concat(
                [
                    qc[qc.test.str.replace('length_', '') == g],
                    qc[qc.test.str.replace('n_percentage_', '') == g],
                    qc[qc.test.str.replace('coverage_', '') == g],
                    qc[qc.test.str.replace('Q1_', '') == g],
                    qc[qc.test.str.replace('median_', '') == g],
                    qc[qc.test.str.replace('Q3_', '') == g]
                ]
            ).reset_index(drop=True)

            if df.size > 0:
                df = df.drop(['ref'], axis=1)
                df = df.T
                df = df[1:]
                df.insert(0, 'protein', g)
                df.columns = ['protein', 'Length', 'N_percentage', 'Coverage', 'Q1_Depth', 'Median_Depth', 'Q3_Depth']
                qc_T.loc[len(qc_T.index)] = df.iloc[0]

        fig2 = make_subplots(
            rows=1, cols=1,
            specs=[[{"type": "table"}]]
        )

        fig2.add_trace(
            go.Table(
                header=dict(values=list(qc_T.columns),
                            fill_color='lightblue',
                            align='left'
                            ),
                cells=dict(values=[
                    qc_T.protein, qc_T.Length, qc_T.N_percentage,
                    qc_T.Coverage, qc_T.Q1_Depth, qc_T.Median_Depth,
                    qc_T.Q3_Depth
                ],
                           fill_color='white',
                           align='left')
                           ),
            row=1, col=1
        )

        fig2.update_layout(title_text="""
            <span style="font-weight: bold;">QC metrics summary</span>
            """
        )
        fig2.write_html(outd_p + '/qc_metrics_summary.html')

    else:
        qc = pd.read_csv(csv, sep=';')
        qc.columns = ["ref", "test", "Count"]

        qc_T = pd.DataFrame(columns=['protein', 'Length', 'Coverage'])

        for g in proteins_list:
            df = pd.concat(
                [
                    qc[qc.test.str.replace('length_', '') == g],
                    qc[qc.test.str.replace('coverage_', '') == g]
                ]
            ).reset_index(drop=True)

            if df.size > 0:
                df = df.drop(['ref'], axis=1)
                df = df.T
                df = df[1:]
                df.insert(0, 'protein', g)
                print(df)
                df.columns = ['protein', 'Length', 'Coverage']
                qc_T.loc[len(qc_T.index)] = df.iloc[0]

        fig2 = make_subplots(
            rows=1, cols=1,
            specs=[[{"type": "table"}]]
        )

        fig2.add_trace(
            go.Table(
                header=dict(values=list(qc_T.columns),
                            fill_color='lightblue',
                            align='left'
                            ),
                cells=dict(values=[qc_T.protein, qc_T.Length, qc_T.Coverage],
                           fill_color='white',
                           align='left')
                           ),
            row=1, col=1
        )

        fig2.update_layout(title_text="""
            <span style="font-weight: bold;">QC metrics summary</span>
            """
        )
        fig2.write_html(outd_p + '/qc_metrics_summary.html')

##### ########################################################################

    muts = csv_muts
    muts.columns = [
        "SampleID", "Protein", "Mutation_type", "Aa_change",
        "Amino_Acid_Property_Change", "Nt_mutation", "Mutation_frequency", "Mutation_depth"
    ]
    muts['ColorCode'] = muts['Mutation_type']
    muts.loc[muts['ColorCode'] == 'NON_SYNONYMOUS', 'ColorCode'] =  '#2243f5'
    muts.loc[muts['ColorCode'] == 'SYNONYMOUS', 'ColorCode'] = '#b9c3fa'
    muts['Type'] = muts['Mutation_type']
    muts.loc[muts['Type'] == 'NON_SYNONYMOUS', 'Type'] =  'Non Synonymous Mutation'
    muts.loc[muts['Type'] == 'SYNONYMOUS', 'Type'] = 'Synonymous Mutation'
    muts['Order'] = muts.Aa_change.str.extract(r'(\d+)', expand=False)
    muts['Order'] = pd.to_numeric(muts['Order'])
    muts = muts.sort_values('Order')

    if not annotate.empty:

        figR = make_subplots(
            rows=rows, cols=1,
            subplot_titles=list(proteins_list)
        )

        muts['Annotated_mutation'] = muts.Aa_change
        muts['Annotated'] = 'no'

        annot_m = annotate
        annot_m.columns = [
            "SampleID", "Protein", "Mutation_type", "Aa_change",
            "Amino_Acid_Property_Change", "Nt_mutation", "Mutation_frequency",
            "Mutation_depth", "Annotated", "Annotated_mutation"
        ]
        annot_m['Order'] = annot_m.Aa_change.str.extract(r'(\d+)', expand=False)
        annot_m['Order'] = pd.to_numeric(annot_m['Order'])
        annot_m = annot_m.sort_values('Order')

        c = 1

        for i in annot_m.Aa_change:
            muts.loc[muts.Aa_change == i, 'Annotated'] = 'yes'
            muts.loc[muts.Aa_change == i, 'Annotated_mutation'] = annot_m.loc[annot_m.Aa_change == i, 'Annotated_mutation'].values[0]
            muts.loc[muts.Aa_change == i, 'ColorCode'] = '#fe0000'
            muts.loc[muts.Aa_change == i, 'Type'] = 'Annotated Mutation'

        muts['TextData'] = muts["Type"] + " = " + muts['Annotated_mutation']
        muts['CustomData'] = "Freq = " + round(muts["Mutation_frequency"], 3).astype(str) + ", " + "Depth = " + muts["Mutation_depth"].astype(str)

        for g in proteins_list:

            muts_p = muts[muts.Protein == g].reset_index(drop=True)

            figR.add_trace(
                go.Scatter(
                    x=muts_p.Order, y=muts_p.Mutation_frequency,
                    mode='markers',
                    name=str(g),
                    text=muts_p.TextData,
                    customdata=muts_p.CustomData,
                    marker_color=muts_p.ColorCode
                ),
                row=c, col=1
            )
            c += 1

        hover_template = "%{text}<br>%{customdata}"
        figR.update_traces(hovertemplate=hover_template)
        for i in range(2, rows + 1):
            figR.update_yaxes(range=[0, 1.1], row=i, col=1)
        for i in range(2, rows + 1):
            end = int(lengths_df.end[lengths_df.row == i].iloc[0])
            figR.update_xaxes(range=[0, end], row=i, col=1)
        figR.update_layout(title_text="""
            <span style="font-weight: bold;">Mutation summary</span>
            """,
            height=400 * rows, showlegend=False
        )
        figR.write_html(outd_p + '/Mutations_summary.html')

        muts_df_end = muts.drop(['ColorCode', 'TextData', 'CustomData', 'Order', 'Type'], axis=1)
        muts_df_end.to_csv(outd_csv + '/' + SAMPLE + '_mutations.csv', sep=';', index=False)
    else:
        fig = make_subplots(
            rows=rows, cols=1,
            subplot_titles=list(proteins_list)
        )

        muts['TextData'] = muts["Type"] + " = " + muts['Aa_change']
        muts['CustomData'] = "Freq = " + round(muts["Mutation_frequency"], 3).astype(str) + ", " + "Depth = " + muts["Mutation_depth"].astype(str)
        c = 1

        for g in proteins_list:

            muts_p = muts[muts.Protein == g].reset_index(drop=True)

            fig.add_trace(
                go.Scatter(
                    x=muts_p.Order, y=muts_p.Mutation_frequency,
                    mode='markers',
                    name=str(g),
                    text=muts_p.TextData,
                    customdata=muts_p.CustomData,
                    marker=dict(color=muts_p.ColorCode)
                ),
                row=c, col=1
            )
            c += 1

        hover_template = "%{text}<br>%{customdata}"
        fig.update_traces(hovertemplate=hover_template)
        for i in range(1, rows + 1):
            fig.update_yaxes(range=[0, 1.1], row=i, col=1)
        for i in range(1, rows + 1):
            end = int(lengths_df.end[lengths_df.row == i].iloc[0])
            fig.update_xaxes(range=[0, end], row=i, col=1)
        fig.update_layout(title_text="""
            <span style="font-weight: bold;">Mutation summary</span>
            """,
            height=400 * (rows + 1), showlegend=False
        )

        fig.write_html(outd_p + '/Mutations_summary.html')
        muts_df_end = muts.drop(['ColorCode', 'TextData', 'CustomData', 'Order', 'Type'], axis=1)
        muts_df_end.to_csv(outd_csv + '/' + SAMPLE + '_mutations.csv', sep=';', index=False)
