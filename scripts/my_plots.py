import os
import seaborn as sns; sns.set(color_codes=True)
import pandas as pd
# import matplotlib.pyplot as plt
from matplotlib import cm
import plotly.express as px
import plotly.io as pio
# from skbio.stats.ordination import pcoa

cm_paired = cm.get_cmap('Set1')
c1 = list(cm_paired.colors)
cm_paired = cm.get_cmap('tab20')
c2 = list(cm_paired.colors)
cm_paired = cm.get_cmap('Set3')
c3 = list(cm_paired.colors)

c_colors = c1 + c3

def crate_graph_folder(output_path):
    # output_path = output_path + "/report"
    output_path = output_path + "/Graphs"
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    print(output_path)
    return output_path


# def create_pie_charts(deneme, output, sample_name):
#     output_path = "/".join(output.split("/")[:-1])
#     output_path = crate_graph_folder(output_path)
#     sample_name_abnd = f"%{sample_name}"
#     levels = list(set(deneme["Taxonomy"].values))
    
#     levels = [level for level in levels if level != "Domain"]

#     level_dict = {} 
#     for level in levels:
#         df = deneme.query("Taxonomy == '{}'".format(level))
#         df.loc[df[sample_name_abnd] < 1, "taxon"] = "Others"
#         df = df.groupby(["Taxonomy", "taxon"], as_index=False).sum()
        
#         pie_df = df.groupby("taxon").sum().reset_index()
#         fig = px.pie(pie_df, names="taxon", values=sample_name_abnd, title=f"{level} Level Taxonomic Comparison of Sample <b>{sample_name}</b>")
#         fig.update_layout(font=dict(size=20), legend_orientation="v", width=1000, height=600)  
        
#         output_path_png = output_path + f"/{sample_name}_piecharts_{level}.png"
#         output_path_html = output_path + f"/{sample_name}_piecharts_{level}.html"
#         fig.write_image(output_path_png)
#         pio.write_html(fig, file=output_path_html, auto_open=False)
        
#         level_dict[level] = output_path_png
#     return level_dict


def create_pie_charts(deneme, output, sample_name):
    output_path = "/".join(output.split("/")[:-1])
    output_path = crate_graph_folder(output_path)
    sample_name_abnd = f"%{sample_name}"
    levels = list(set(deneme["Taxonomy"].values))
    
    levels = [level for level in levels if level != "Domain"]

    level_dict = {} 
    for level in levels:
        df = deneme.query("Taxonomy == '{}'".format(level))
        df.loc[df[sample_name_abnd] < 1, "taxon"] = "Others"
        df = df.groupby(["Taxonomy", "taxon"], as_index=False).sum()
        
        pie_df = df.groupby("taxon").sum().reset_index()
        fig = px.pie(pie_df, names="taxon", values=sample_name_abnd, title=f"{level} Level Taxonomic Comparison of Sample <b>{sample_name}</b>")
        
        # Güncellenmiş layout ayarları
        fig.update_layout(
            font=dict(size=20),
            legend_orientation="v", 
            width=1000, 
            height=600,
            legend=dict(
                x=1.1,  # x ekseninde legend'in pozisyonu (chart'tan sağa kaydırır)
                y=0.5,  # y ekseninde legend'in ortalanması
                xanchor="left"  # legend'in hizalaması
            )
        )
        
        output_path_png = output_path + f"/{sample_name}_piecharts_{level}.png"
        output_path_html = output_path + f"/{sample_name}_piecharts_{level}.html"
        fig.write_image(output_path_png)
        pio.write_html(fig, file=output_path_html, auto_open=False)
        
        level_dict[level] = output_path_png
    return level_dict



def create_per_seq_quality_chart(df, output, sample_name):
    output_path = "/".join(output.split("/")[:-1])
    output_path = crate_graph_folder(output_path)
    
    df = pd.DataFrame(df, columns=['Quality', 'Count'], dtype=float)
    
    fig = px.line(df, x='Quality', y='Count', title=f'Quality vs. Count for Sample <b>{sample_name}</b>')
    fig.update_traces(mode='markers+lines')
    fig.update_layout(font=dict(size=20), plot_bgcolor='white', width=1000, height=600, 
                      xaxis=dict(
                                tickmode='array',
                                tickvals=df['Quality'].tolist(), 
                                ticktext=[str(int(val)) for val in df['Quality'].tolist()],
                                title='Mean Sequence Quality (Phred Score)',
                                showgrid=True,
                                showline=True,
                                zeroline=True,
                                zerolinecolor='#bdbdbd',
                                zerolinewidth=2,
                                mirror='ticks',
                                gridcolor='#bdbdbd',
                                gridwidth=2
                    ),
                    yaxis=dict(
                                title='Count',
                                showgrid=True,
                                showline=True,
                                mirror='ticks',
                                zeroline=True,
                                zerolinecolor='#bdbdbd',
                                zerolinewidth=2,
                                gridcolor='#bdbdbd',
                                gridwidth=2
                    ),
                    margin=dict(l=50, r=50, t=50, b=50),
                    shapes=[
                            dict(
                                type="rect",
                                xref="paper",
                                yref="paper",
                                x0=0,
                                y0=0,
                                x1=1,
                                y1=1,
                                line=dict(
                                    color="black",
                                    width=2
                                )
                            )
                        ])
    fig.update_traces(marker_size=10)
    sub_output = output_path + f"/{sample_name}_linegraph.png"    
    output_path_ = output_path + f"/{sample_name}_linegraph.html"

    fig.write_image(sub_output)
    pio.write_html(fig, file=output_path_, auto_open=False)
    return sub_output, output_path_


def quality_plot(output, df):
    o_file_name = output.split("/")[-1].split(".")[0]
    output_path = "/".join(output.split("/")[:-2])
    output_path = output_path + "/report"
    output_path = crate_graph_folder(output_path)
    print(df)
    ax = sns.jointplot(df["normal_avg"], df["read_lenght"], kind="kde", height=7, space=0,xlim=(15,25), ylim=(1300,1800))
    ax.ax_joint.set_xlabel('Quality Scores', fontsize = 15)
    ax.ax_joint.set_ylabel('Read Length', fontsize = 15)
    output_path = output_path + "/quality_plot_{}.png".format(o_file_name)
    ax.savefig(output_path)
    return output_path
