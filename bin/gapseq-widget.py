#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import plotly.express as px
from glob import glob
from scipy.cluster import hierarchy
import widgets.streamlit as wist


class GapseqWidget(wist.StreamlitWidget):

    title = 'Microbial Genome Metabolism Explorer'
    
    disable_sidebar = True

    extra_imports = [
        "import plotly.express as px",
        "from glob import glob",
        "from scipy.cluster import hierarchy"
    ]

    requirements = ["scipy", "plotly"]

    children = [
        wist.StExpander(
            id='options',
            children=[
                wist.StDataFrame(
                    id='pathways',
                    show_uploader=False
                ),
                wist.StMultiSelect(
                    id='show_samples',
                    label='Show Samples',
                    options=[]
                ),
                wist.StInteger(
                    id='height',
                    min_value=200,
                    max_value=2500,
                    value=700
                ),
                wist.StInteger(
                    id='width',
                    min_value=200,
                    max_value=2500,
                    value=700
                )
            ]
        )
    ]

    def read_inputs(self, basedir='.', suffix='-all-Pathways.tbl'):

        all_inputs = []
        all_samples = []

        for fp in glob(f'{basedir}/**/*{suffix}', recursive=True):

            # Read in the pathways
            df = pd.read_csv(
                fp,
                sep='\t',
                skiprows=3
            )

            # Filter down to just the predicted ones
            df = df.query('Prediction')

            # Parse the sample name
            sample_name = fp.rsplit('/', 1)[-1][:-len(suffix)]
            all_samples.append(sample_name)

            # Add it to the table
            df = df.assign(
                sample_name=sample_name,
                ID=df['ID'].apply(lambda s: s.strip('|'))
            )

            all_inputs.append(df)

        all_inputs = pd.concat(all_inputs)

        self.set(path=['options', 'pathways'], value=all_inputs, update=False)
        self.set(path=['options', 'show_samples'], attr='options', value=all_samples, update=False)
        self.set(path=['options', 'show_samples'], attr='value', value=all_samples, update=False)

    def run_self(self, **kwargs) -> None:

        # Get the set of pathways
        df_long = self.get(['options', 'pathways'])

        # Filter to the set of samples which were indicated
        df_long = df_long.loc[
            df_long['sample_name'].isin(self.get(['options', 'show_samples']))
        ]

        # Make a wide presence-absence table
        df_wide = df_long.pivot_table(
            index='sample_name',
            columns='Name',
            values='Completeness'
        ).fillna(0)

        # Order by linkage
        df_wide = df_wide.reindex(
            index=self.order_index_by_linkage(df_wide),
            columns=self.order_index_by_linkage(df_wide.T)
        )
        
        # Make a heatmap
        fig = px.imshow(
            df_wide,
            aspect='auto',
            height=self.get(['options', 'height']),
            width=self.get(['options', 'width']),
            labels=dict(
                color='Completeness (%)'
            )
        )

        fig.update_layout(
            xaxis_title='',
            yaxis_title=''
        )

        self.main_container.plotly_chart(fig)

        self.download_html_button()
        self.download_script_button()

    def order_index_by_linkage(self, df):

        return df.index.values[
            hierarchy.leaves_list(
                hierarchy.linkage(
                    df.values
                )
            )
        ]


if __name__ == '__main__':

    w = GapseqWidget()

    w.read_inputs()

    # w.run()
    w.to_html(Path("gapseq_pathways.html"))
