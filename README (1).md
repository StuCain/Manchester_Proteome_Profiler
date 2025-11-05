
# Manchester Proteome Profiler (MPP)

Manchester Proteome Profiler (MPP) is an open-source R Shiny application designed for downstream analysis of quantitative proteomics data. It supports input from MaxQuant, FragPipe, Proteome Discoverer, and Spectronaut, and provides a user-friendly interface for data filtering, normalisation, imputation, differential expression analysis, clustering, and functional enrichment.

---

## 🚀 Features
- Supports single and dual dataset comparisons
- Normalisation and imputation options
- Differential expression analysis using limma
- SAINTexpress support for AP-MS/BioID data
- Gene Ontology enrichment and protein interaction networks (BioGRID, STRING)
- Interactive visualisations: PCA, volcano plots, heatmaps, violin plots
- Publication-ready plots and downloadable reports

---

## 🐳 Docker Deployment
To run MPP using Docker:

```bash
docker pull stu-cain/mpp
docker run -d -p 3838:3838 stu-cain/mpp
```

Then open your browser and go to: [http://localhost:3838](http://localhost:3838)

---

## 📦 Installation (Manual)
1. Clone the repository:
   ```bash
   git clone https://github.com/StuCain/Manchester_Proteome_Profiler.git
   ```
2. Install R and required packages (see Dockerfile for dependencies)
3. Launch the app:
   ```R
   shiny::runApp('path/to/shiny-mpp')
   ```

---

## 📁 Input Requirements
- Grouped protein intensity files (e.g., `proteinGroups.txt`, `combined_protein.tsv`)
- Experimental design file (template provided in-app)
- Optional: comparison file for complex queries

---

## 📚 Documentation
Full documentation and usage examples are available in the GitHub Wiki and within the app interface.

---

## 📖 Citation
If you use MPP in your research, please cite:
> Cain, S.A. et al. Manchester Proteome Profiler: A User-Friendly Platform for Quantitative Proteomics Analysis. *Journal of Proteome Research*, 2025.

---

## 🛠 License
This project is licensed under the MIT License.
