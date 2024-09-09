#!/usr/bin/env python3
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np


def visualizer(tsv_check, csv, csv_muts, outd_p, outd_csv, genes, annotate, SAMPLE, syn_muts, lengths_df):
    genes_list = genes
    rows = len(genes_list)

    if tsv_check:
        qc = pd.read_csv(csv, sep=';')
        qc.columns = ["ref", "test", "Count"]

        qc_T = pd.DataFrame(columns=['gene', 'Length', 'N_percentage', 'Coverage', 'Q1_Depth', 'Median_Depth', 'Q3_Depth'])

        for g in genes_list:
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
                df.insert(0, 'gene', g)
                df.columns = ['gene', 'Length', 'N_percentage', 'Coverage', 'Q1_Depth', 'Median_Depth', 'Q3_Depth']
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
                    qc_T.gene, qc_T.Length, qc_T.N_percentage,
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

        qc_T = pd.DataFrame(columns=['gene', 'Length', 'Coverage'])

        genes_list = genes

        rows = len(genes_list)

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
                cells=dict(values=[qc_T.gene, qc_T.Length, qc_T.Coverage],
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
        "SampleID", "Gene", "Mutation_type", "Aa_change",
        "Type_of_aa_change", "Nt_mutation", "Mutation_frequency"
    ]
    muts['ColorCode'] = muts['Mutation_type']
    muts.loc[muts['ColorCode'] == 'MUTATION', 'ColorCode'] =  '#2243f5'
    muts.loc[muts['ColorCode'] == 'NO_MUTATION', 'ColorCode'] = '#b9c3fa'
    muts['Type'] = muts['Mutation_type']
    muts.loc[muts['Type'] == 'MUTATION', 'Type'] =  'Non Synonymous Mutation'
    muts.loc[muts['Type'] == 'NO_MUTATION', 'Type'] = 'Synonymous Mutation'
    muts['Order'] = muts.Aa_change.str.extract(r'(\d+)', expand=False)
    muts['Order'] = pd.to_numeric(muts['Order'])
    muts = muts.sort_values('Order')

    if not annotate.empty:

        figR = make_subplots(
            rows=rows, cols=1,
            subplot_titles=list(genes_list)
        )

        muts['Annotated_mutation'] = muts.Aa_change
        muts['Annotated'] = 'no'

        annot_m = annotate
        annot_m.columns = [
            "SampleID", "Gene", "Mutation_type", "Aa_change",
            "Type_of_aa_change", "Nt_mutation", "Mutation_frequency",
            "Annotated", "Annotated_mutation"
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
        muts['CustomData'] = "Freq = " + round(muts["Mutation_frequency"], 3).astype(str)

        for g in genes_list:

            muts_p = muts[muts.Gene == g].reset_index(drop=True)

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
        # Customize axes for each subplot
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

        muts_df_end = muts.drop(['ColorCode', 'TextData', 'CustomData', 'Order'], axis=1)
        muts_df_end.to_csv(outd_csv + '/' + SAMPLE + '_mutations.csv', sep=';', index=False)
    else:
        fig = make_subplots(
            rows=rows, cols=1,
            subplot_titles=list(genes_list)
        )

        muts['TextData'] = muts["Type"] + " = " + muts['Aa_change']
        muts['CustomData'] = "Freq = " + round(muts["Mutation_frequency"], 3).astype(str)
        c = 1

        for g in genes_list:

            muts_p = muts[muts.Gene == g].reset_index(drop=True)

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
        # Customize axes for each subplot
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
