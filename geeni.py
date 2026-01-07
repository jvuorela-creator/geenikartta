import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import numpy as np

# Asetukset ja sivun konfiguraatio
st.set_page_config(page_title="Genomi-kartografia", layout="wide", page_icon="🧬")

def load_data(uploaded_file):
    """Lataa FTDNA Chromosome Browser -tiedoston."""
    try:
        # FTDNA Chromosome Browser tiedostossa on yleensä otsikot:
        # "First Name", "Last Name", "Match Name", "Chromosome", "Start Location", "End Location", "Centimorgans", "Matching SNPs"
        df = pd.read_csv(uploaded_file)
        
        # Nimien siistiminen (poistetaan välilyönnit sarakkeiden nimistä)
        df.columns = df.columns.str.strip()
        
        # Varmistetaan, että tarvittavat sarakkeet löytyvät
        required_cols = ['Chromosome', 'Start Location', 'End Location', 'Centimorgans', 'Match Name']
        if not all(col in df.columns for col in required_cols):
            st.error("Tiedostosta puuttuu vaadittuja sarakkeita. Lataa FTDNA Chromosome Browser -CSV.")
            return None
            
        return df
    except Exception as e:
        st.error(f"Virhe tiedoston luvussa: {e}")
        return None

def create_topography(df):
    """Luo 3D-maaston DNA-segmenteistä."""
    
    # Alustetaan kuvaaja
    fig = go.Figure()

    # Määritellään väripaletti maastolle (Meri -> Ranta -> Metsä -> Vuori -> Lumi)
    colors = ['#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # Käydään läpi kromosomit 1-22 (X-kromosomi voidaan lisätä haluttaessa)
    chromosomes = [str(i) for i in range(1, 23)]
    
    # Etsitään uniikit osumat (serkut)
    matches = df['Match Name'].unique()
    
    # Jos osumia on paljon, otetaan vain top 5 selkeyden vuoksi, ellei käyttäjä valitse toisin
    if len(matches) > 10:
        st.warning(f"Tiedostossa on {len(matches)} osumaa. Näytetään vain 5 suurinta selkeyden vuoksi.")
        # Lasketaan kokonais-cM jokaiselle ja otetaan top 5
        match_totals = df.groupby('Match Name')['Centimorgans'].sum().sort_values(ascending=False)
        matches = match_totals.head(5).index.tolist()
        df = df[df['Match Name'].isin(matches)]

    max_position = 250_000_000 # Arvioitu maksimipituus kromosomille visualisointia varten
    
    # Luodaan "harjanne" jokaiselle kromosomille
    for i, chrom in enumerate(chromosomes):
        chrom_data = df[df['Chromosome'] == chrom]
        
        # Luodaan perusviiva (merenpinta)
        x_base = np.linspace(0, max_position, 500)
        y_base = np.full_like(x_base, 0)
        
        # Jos kromosomissa on dataa, muokataan korkeutta
        if not chrom_data.empty:
            for _, row in chrom_data.iterrows():
                start = row['Start Location']
                end = row['End Location']
                cm = row['Centimorgans']
                match_name = row['Match Name']
                
                # Luodaan "vuori" segmentin kohdalle
                # Käytetään maskia valitsemaan oikeat kohdat x-akselilla
                mask = (x_base >= start) & (x_base <= end)
                
                # Lisätään korkeutta (cM) niihin kohtiin
                # Pieni satunnaisuus tai kerrostaminen auttaa erottamaan päällekkäiset osumat
                y_base[mask] += cm 

        # Lisätään täytetty alue (vuoristo)
        # Y-akseli on tässä tapauksessa kromosomin numero (väännettynä 3D-tilaan)
        # Z-akseli on korkeus (cM)
        
        # Koska Plotly 3D vaatii hieman kikkailua "Ridge Plot" tyyliin:
        # Teemme jokaisesta kromosomista oman "viivan" Y-akselilla
        
        fig.add_trace(go.Scatter3d(
            x=x_base,
            y=np.full_like(x_base, 23 - int(chrom)), # Käännetään järjestys: Chr 1 ylhäällä
            z=y_base,
            mode='lines',
            line=dict(color='black', width=1),
            name=f'Chr {chrom}',
            showlegend=False
        ))

        # Pintaefekti (täytetään viivan alapuoli)
        # Tämä on visuaalinen temppu: piirretään pinta alas nollaan
        # 3D:ssä tämä on hieman raskasta, joten käytämme viivoja "topografisena karttana"
        # Vaihtoehto: Surface plot interpoloimalla, mutta viivat ovat tarkempia geneettisesti.
        
        # Värjätään korkeuden mukaan
        color_val = y_base.max()
        if color_val > 20: c = 'firebrick' # Iso vuori
        elif color_val > 10: c = 'forestgreen' # Kukkula
        elif color_val > 0: c = 'sandybrown' # Saari
        else: c = 'aliceblue' # Meri
        
        # Piirretään "täyte" pystysuorilla viivoilla (kuten aita) tietyin välein
        # jotta se näyttää kiinteältä maastolta
        subset_idx = np.where(y_base > 0)[0]
        if len(subset_idx) > 0:
             fig.add_trace(go.Scatter3d(
                x=x_base[subset_idx],
                y=np.full_like(x_base[subset_idx], 23 - int(chrom)),
                z=y_base[subset_idx],
                mode='markers',
                marker=dict(size=2, color=c),
                showlegend=False,
                hoverinfo='text',
                text=[f"Chr {chrom}: {val:.1f} cM" for val in y_base[subset_idx]]
            ))

    fig.update_layout(
        title="Genomi-Topografia: Serkkuosumat maisemana",
        scene=dict(
            xaxis_title='Sijainti kromosomissa (bp)',
            yaxis_title='Kromosomi (1-22)',
            zaxis_title='Yhteyden voimakkuus (cM)',
            yaxis=dict(tickvals=list(range(1, 23)), ticktext=[f"Chr {23-i}" for i in range(1, 23)]),
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=0.5) # Viisto ilmakuva
            )
        ),
        margin=dict(l=0, r=0, b=0, t=50)
    )
    
    return fig

# --- KÄYTTÖLIITTYMÄ ---

st.title("🧬 Genomi-kartografia: DNA-Topografia")

st.markdown("""
**Tervetuloa sukututkija.** Tämä työkalu muuttaa FamilyTreeDNA:n datan maisemaksi.
* **Vuoristot:** Suuret yhteiset segmentit (läheinen sukulaisuus).
* **Saaret:** Pienet, satunnaiset tai kaukaiset osumat.
* **Laaksot:** Alueet, joissa ei ole yhteistä perimää.
""")

col1, col2 = st.columns([1, 2])

with col1:
    st.info("""
    **Ohje:** Lataa tähän **Chromosome Browser Results** -tiedosto (CSV).
    
    Älä lataa "Build 37 Raw Data" -tiedostoa, sillä se sisältää vain sinun perimäsi kirjaimet. 
    Maiseman luomiseksi tarvitsemme *segmenttejä*, jotka syntyvät vertailusta osumiin.
    """)
    uploaded_file = st.file_uploader("Lataa CSV-tiedosto", type=['csv'])

if uploaded_file is not None:
    df = load_data(uploaded_file)
    
    if df is not None:
        with col2:
            st.success(f"Ladattu {len(df)} segmenttiä.")
            
            # Suodattimet
            min_cm = st.slider("Minimi cM (suodata kohinaa)", 1, 20, 5)
            df_filtered = df[df['Centimorgans'] >= min_cm]
            
            st.write(f"Näytetään segmentit, jotka ovat suurempia kuin {min_cm} cM.")
        
        # Visualisointi
        st.subheader("DNA-Maisema")
        fig = create_topography(df_filtered)
        st.plotly_chart(fig, use_container_width=True)
        
        # Datataulukko
        with st.expander("Katso lähdedata"):
            st.dataframe(df_filtered)

else:
    # Demo-näkymä jos dataa ei ole
    st.markdown("---")
    st.write("*Lataa tiedosto nähdäksesi visualisoinnin.*")
