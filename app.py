from shiny import App, ui, render, reactive, run_app
from shiny.ui import update_text_area
import threading
import time
import webbrowser
import pandas as pd
from htmltools import tags
from evercpt import main, table_motifs, plot_motifs

# UI
app_ui = ui.page_fluid(
    ui.br(),
    tags.head(
        tags.style(
            """
            body {
                background-color: #ffffff;
                color: #000000;
                transition: background-color 0.5s, color 0.5s;
                margin: 0;
                padd()ing: 0;
            }
            .btn-custom {
                background-color: #007bff;
                color: white;
                border: none;
                border-radius: 4px;
                padd()ing: 8px 16px;
                font-size: 14px;
                display: block;
                margin: 0 auto;
            }
            .btn-custom2 {
                background-color: #999999;
                color: white;
                border: none;
                border-radius: 4px;
                padding: 8px 16px;
                font-size: 14px;
                display: block;
                width: 100%;
            }
            .btn-custom:hover {
                background-color: #0056b3;
            }
            .text-input {
                height: 100px;
                width: 100%;
                box-sizing: border-box;
            }
            table {
                border-collapse: collapse;
            }
            table tr:nth-child(odd) {
                background-color: #f2f2f2;
            }
            table.shiny-table tr:nth-child(odd) {
                background-color: #f2f2f2 !important;
            }
            .btn-custom:hover {
                background-color: #0056b3;
            }
            .text-input {
                height: 100px;
                width: 100%;
                box-sizing: border-box;
            }
            .container-fluid {
                width: 100%;
                max-width: 800px;
                margin: 0 auto;
                padding-top: 20px;
            }
            .bottom-tabs .nav-tabs {
                display: flex;
                justify-content: space-between;
                flex-wrap: wrap;
                width: 100%;            /* full width */
                border-bottom: 1px solid #ccc;  /* line directly under tabs */
                margin: 0;
                padding: 0 2.5%;       /* padding left+right = 25%, centers tabs inside */
                box-sizing: border-box;
            }
            .bottom-tabs .nav-tabs > li {
                flex: 1;
                text-align: center;
            }
            .bottom-tabs .nav-tabs > li > a {
                width: 100%;
            }
            .bottom-tabs .tab-content {
                margin-top: 1rem;
            }
            th, td {
                text-align: left !important;
            }
            """
        )
    ),
    
    ui.h1("EVERCPT: Extracellular Vesicles Enrichment of RNA Cargo Prediction Tool", style="text-align:center;"),
    ui.navset_tab(
        
        ui.nav_panel(
            
            "Custom", ui.br(),
    
            ui.input_select("dropdown", "Select a species:", choices={"human": "Human"}, width="200px"),

            ui.input_radio_buttons("radio", "Choose an option:", choices={"mrna": "mRNA", "circ": "circRNA"}),

            ui.row(
                ui.column(12,  # full width
                    ui.tags.div(
                        ui.input_action_button("example_btn","Example sequence",class_="btn-custom2",style="width: 160px;"
                        ),
                        ui.input_action_button("clear_btn","Clear",class_="btn-custom2",style="width: 100px;"
                        ),
                        style="display: flex; justify-content: flex-end; gap: 12px; margin-top: 30px;"
                    )
                )
            ),

            ui.input_text_area("textInput", "Enter sequence:", value="", width="800px", rows=5),
            
            ui.row(ui.column(12, ui.input_action_button("extractButton", "Predict", class_="btn-custom"))),

            ui.br(),
            ui.output_code("text_result"),
            ui.br(),

            ui.div(
                ui.div(
                    ui.tags.div(class_="tabs-header"),
                    ui.navset_tab(
                        ui.nav_panel("Features",
                            ui.output_table("table1", striped=True, bordered=True, hover=True),
                            ui.output_table("table2", striped=True),
                            ui.output_table("table3", striped=True),
                            ui.output_table("table4", striped=True),
                        ),
                        ui.nav_panel("Motifs",
                                ui.br(),
                                ui.layout_columns(
                                    ui.div(
                                        ui.output_table("motif_table", striped=True),
                                        style="width: 100%; max-width: 500px; max-height: 700px; overflow-x: auto; overflow-y: auto;"
                                    ),
                                    ui.div(
                                        ui.output_plot("motif_plot", width="100%", height="300px"),
                                        style="height: 400px; overflow: hidden;"
                                    ),
                                    col_widths=[6, 6]
                                )
                        ),
                    ),
                    class_="bottom-tabs"
                ),
                class_="bottom-tabs-wrapper"
            )

        ),
        ui.nav_panel(""),

    )
)

# SERVER
def server(input, output, session):
    @reactive.effect
    @reactive.event(input.example_btn)
    def _():
        update_text_area("textInput", value="CAATGATGTTGTCCACTGGGCATGTACTGACCAATGTGGCAGGTCTGAGAACATAGCTGAAGCTGAAAATAGGAAAGCT"
                         + "GGGGGCAAGGAAGAGCCTTGAATCTTGAGGTGGGACGTTGACTCTAAGATGTCCTTGAGCAGTGGAGCCTCCGGAGGGAAAGGAGTGGATGCAAAC"
                         + "CCGGTTGAGACATACGACAGTGGGGATGAATGGGACATTGGAGTAGGGAATCTCATCATTGACCTGGACGCCGATCTGGAAAAGGACCAGCAGAAA"
                         + "CTGGAAATGTCAGGCTCAAAGGAGGTGGGGATACCGGCTCCCAATGCTGTGGCCACACTACCAGACAACATCAAGTTTGTGACCCCAGTGCCAGGT"
                         + "CCTCAAGGGAAGGAAGGCAAATCAAAATCCAAAAGGAGTAAGAGTGGCAAAGACACTAGCAAACCCACTCCAGGGACTTCCCTGTTCACTCCAAGT"
                         + "GAGGGGGCAGCTAGCAAGAAAGAGGTGCAGGGGCGCTCAGGAGATGGTGCCAATGCTGGAGGCCTGGTTGCTGCTATTGCTCCCAAGGGCTCAGAG"
                         + "AAGGCGGCTAAGGCATCCCGCAGTGTAGCCGGTTCCAAAAAGGAGAAGGAGAACAGCTCATCTAAGAGCAAGAAGGAGAGAAGCGAAGGAGTGGGG"
                         + "ACTTGTTCAGAAAAGGATCCTGGGGTCCTCCAGCCAGTTCCCTTGGGAGGACGGGGTGGTCAGTATGATGGAAGTGCAGGGGTGGATACAGGAGCT"
                         + "GTGGAGCCACTTGGGAGTATAGCTATTGAGCCTGGGGCAGCGCTCAATCCTTTGGGAACTAAACCGGAGCCAGAGGAAGGGGAGAATGAGTGTCGC"
                         + "CTGCTAAAGAAAGTCAAGTCTGAAAAG")

    @reactive.effect
    @reactive.event(input.clear_btn)
    def _():
        update_text_area("textInput", value="")

    @reactive.Calc
    @reactive.event(input.extractButton)
    def result():
        seq = input.textInput().strip()
        if seq:
            return main(seq=seq, type=input.radio())

    
    @reactive.Calc
    @reactive.event(input.extractButton)
    def add():
        return 1 if input.radio() == "circ" else 0
    
    @output
    @render.code
    @reactive.event(input.extractButton)
    def text_result():
        res = result()
        if res is not None:
            score = res["result"].iloc[0]
            if score > 0.5:
                return f"This RNA is likely enriched in EVs (enrichment score: {score:.3f})"
            else:
                return f"This RNA is likely retained in cells (enrichment score: {score:.3f})"
        else:
            return "Error: Please enter a sequence."


    @output
    @render.table
    def table1():
        res = result()
        if res is not None:
            res = res.iloc[:, :9+add()].copy()                                         # Display first 7 columns as features
            res.loc[:,'Length'] = res['Length'].astype(int)                     # Convert Length to integer
            res.loc[:,'MFE'] = res['MFE'] * -1
            res.loc[:, ['GC%', 'AT%']] = (res.loc[:, ['GC%', 'AT%']]).round(2)
            return res.round(4)                                               # Round all values to 2 decimal places                    

    @output
    @render.table
    def table2():
        res = result()
        if res is not None:
            return (res.iloc[:, 9+add():13+add()]).round(2)                      # Display nucleotide frequencies
    
    @output
    @render.table
    def table3():
        res = result()
        if res is not None:
            return (res.iloc[:, 13+add():21+add()]).round(2)                  # Display dinucleotide frequencies
    
    @output
    @render.table
    def table4():
        res = result()
        if res is not None:
            return (res.iloc[:, 21+add():29+add()]).round(2)                    # Display dinucleotide frequencies

    @output
    @render.table
    def motif_table():
        seq = input.textInput().strip()
        if seq and len(seq) > 10:
            res = table_motifs(seq=seq)
            return res.rename(columns={'sum_motif': 'Distinct Motifs', 'sum_count': 'Motif Count'}).drop(columns=['sum_all'])

    @output
    @render.plot
    def motif_plot():
        seq = input.textInput().strip()
        if seq and len(seq) > 10:
            return plot_motifs(seq=seq)

app = App(app_ui, server)

def open_browser():
    time.sleep(1)
    webbrowser.open("http://localhost:8000")


# Run the app
if __name__ == "__main__":
    threading.Thread(target=open_browser).start()
    run_app(app, port=8000, host="0.0.0.0")
