

**A Machine Learning Framework for CO₂ Activation by Frustrated Lewis Pairs: Database Curation, Prediction, and Discovery**

**Team Name: Error 404**

**Team Members: Kishore Jaykumar Ishwar, Soham Joshi, Aditya Khatu**

# **Project Summaray:**

## **Overview:**

Rising atmospheric carbon dioxide (CO₂) from human activities is a major driver of global warming. Addressing this challenge requires integrated solutions within the carbon capture, storage, and utilization (CCSU) framework—not only capturing CO₂ but also converting it into value-added products. One promising pathway involves **Frustrated Lewis Pairs (FLPs)**—unique combinations of Lewis acids and bases that, because of steric hindrance, remain “frustrated” yet highly reactive. These metal-free systems can heterolytically split hydrogen (H₂) and activate CO₂ under mild conditions, providing a tunable and sustainable platform for catalytic CO₂ conversion..

A recent breakthrough by **Ye et al. (2025)** introduced the first large-scale FLP Database (FLPDB), which contains more than 100,000 acid–base pairs with quantum-chemical descriptors and computed activation data for both H₂ and CO₂. Building on this foundation, the present project focuses specifically on the CO₂ activation subset of that database. We will curate and refine CO₂-specific FLP entries, train machine learning (ML) models to predict how different acid–base combinations interact with CO₂, and apply active-learning strategies to discover new, more efficient FLPs.

By combining modern data science with molecular chemistry, this project aims to accelerate the identification of next-generation, metal-free catalysts for CO₂ utilization. Expected outcomes include a curated, openly available FLP–CO₂ dataset; benchmark ML models ranging from interpretable baselines to advanced graph neural networks; and a shortlist of promising new FLP candidates predicted to capture and convert CO₂ efficiently. 

## **Intellectual Merit:**

This project sits at the intersection of computational chemistry, cheminformatics, and machine learning. It leverages the Ye et al. FLP database to perform the first focused ML analysis on CO₂ activation. Through systematic curation, feature engineering, and modeling, we will build predictive frameworks that reveal how molecular features such as Lewis acidity, base strength, and steric effects influence CO₂ activation.

The work advances the emerging field of molecular machine learning for catalysis by demonstrating how open databases and AI-driven methods can generate reliable chemical predictions. Beyond producing new candidate FLPs, the project will yield broadly applicable insights and computational workflows for other small-molecule activation reactions. We will gain hands-on experience in data processing, molecular modeling, and AI methods, fostering interdisciplinary expertise.

**Broader Impacts:**

This project advances both sustainability and education. Scientifically, it contributes to the global effort to develop metal-free catalysts for CO₂ capture and conversion—an essential step toward a circular carbon economy. By accelerating catalyst discovery with machine learning, it supports technologies that can reduce carbon emissions and reliance on fossil fuels.

Educationally, the project trains students in data-driven chemistry, bridging machine learning and molecular science. All datasets, code, and models will be released openly to promote transparency, reproducibility, and community engagement. The resulting FLP–CO₂ database will serve as both a research resource and a teaching tool, fostering collaboration and innovation in sustainable chemical design.

# **Project Description:** 

## **Problem Framing**

The global Grand Challenge is the escalating atmospheric accumulation of CO2​, necessitating the development of efficient, cost-effective catalytic strategies for its conversion into value-added chemicals.

Frustrated Lewis Pairs (FLPs)—metal-free systems that activate small molecules through cooperative acid-base interactions—represent a promising, inexpensive alternative to traditional precious metal catalysts for CO2​ utilization. However, the rational design of optimal FLPs is severely hindered by the vast, combinatorial chemical space of Lewis acid and base combinations. Systematic structure-property relationships for FLP–CO2​ binding are fragmented, and an integrated, standardized, and predictive computational framework has been lacking, making exhaustive experimental screening impractical.

The recent publication of the **FLPDB (Ye et al., 2025\)** provides the critical opportunity to overcome this design bottleneck. The database aggregates computed properties for over 100k FLPs, including DFT-optimized CO2​ binding energies for a benchmark set of structures. This resource finally makes systematic Machine Learning (ML) analysis feasible.

## **Project Objectives**

The core problem to be solved is the lack of a reliable, high-throughput method to predict and rank the CO2​ binding affinity of novel FLPs. This project addresses this by defining two integrated objectives:

1. Objective 1: Dataset Construction and Standardization  
   * Build a CO2​-focused, open dataset of FLPs by ensuring consistent thermochemistry across the entire corpus. All free energies (ΔG) and relevant electronic/steric descriptors must be computed and standardized under a single Density Functional Theory (DFT) setting. This rigor is essential to ensure that all energy and feature values are comparable for ML training.  
2. Objective 2: Predictive Modeling and Ranking  
   * Test and train ML models (specifically Ridge, Lasso, and Kernel Ridge regression) to accurately predict the CO2​ binding energy (ΔG) of novel FLPs.  
   * Rank new FLPs based on these predictive metrics and the associated uncertainty of the prediction.

## **Datasets / Provenance**

This project will create a new, proprietary dataset for which the provenance will be entirely computational. Our methodology is directly inspired by the high-throughput Density Functional Theory (DFT) screening approach used by Ye et al. to create the Frustrated Lewis Pairs Database (FLPDB). Since the FLPDB is not accessible, we will generate our own dataset by collecting FLP structures from various relevant research papers. This new database will be strategically designed for our research goal of CO2​ activation. All data points will be produced using consistent and validated DFT protocols to ensure high quality and reproducibility. The dataset will be populated with key performance metrics, including calculated binding free energies and activation energy barriers, forming a robust foundation for training our predictive machine learning models.

## **Baselines**

**Model Baselines** \- Train and evaluate Ridge, Lasso, and GBDT on the full feature set for each FLP (geometry, sterics, electronics). Report RMSE, MAE, and R², and assess uncertainty via model ensembles (gives a better estimate of Gibbs free binding energy).

## **Method**

1) Curation and Geometry-  
   We first assemble a set of FLPs by pairing Lewis acids (LAs) and Lewis base (LBs), then build a 3D structure for unbounded FLPs, CO₂, and the FLP.CO₂ adduct. DFT optimization is a local, gradient \- based search on the potential energy surface. To avoid getting stuck in sub-optimal structure, try multiple placements of CO₂ relative to the lewis acid and base, optimize each case, then conclude the lowest-energy result. DFT—greatly reducing the chance of converging to a poor local minimum that would distort binding energies.  
     
2) Thermochemistry-  
   Using one DFT setting (functional,basis,solvation model) for every system, we optimize the geometries in the gas phase and run frequency analysis to obtain thermodynamic properties like enthalpy, entropy which yields gas \- phase Gibbs free energy. 

3) Label-  
   Label indicates how strongly does CO₂ binds with a specific FLP in solution and label is the factor which the model needs to learn to predict from inputs. During training, the model takes input which is features of each FLP pair—geometry, steric, and electronic descriptors—along with the label, the solution-phase Gibbs free energy of binding. It learns a mapping from features to Gibbs free energy of binding that generalizes to new FLPs. At inference, the model predicts Gibbs free energy of binding for unseen candidates; because more negative \= stronger binding, we rank FLPs by the predicted label.  
     
4) Metrics-  
   For model evaluation we report the regression metrics:  
* RMSE: overall prediction error magnitude.  
* MAE: average absolute error (robust to outliers).  
* R²: fraction of variance explained.


  For decision-making we rank candidates by predicted ΔGbind (more negative is better) and select the top entries for confirmatory DFT.

### **6\. Risks and Mitigation**

| Risk | Potential Impact | Mitigation Strategy |
| :---- | ----- | ----- |
| **Limited access to FLPDB or insufficient data points** | Could delay dataset generation | Use published FLP studies as alternative sources; automate DFT data generation through in-house scripts; prioritize a smaller, high-quality dataset. |
| **DFT data generation is computationally expensive** | Delayed data generation and model training | build a **predictive model** to approximate expensive computations, then **iteratively select and validate** only the most informative or uncertain data points, reducing total simulation cost while improving model accuracy. |
| **Model overfitting due to small dataset** | Reduces predictive reliability | Use **simpler models with regularization**, apply **cross-validation** to check consistency, and **expand data** through active learning. |
| **Label noise and convergence errors** | Introduces data inaccuracies | Exclude unstable geometries (imaginary frequencies \> 1); repeat optimization; validate subset against higher-level DLPNO-CCSD(T) single points. |
| **Time management** | Risk of incomplete deliverables | Weekly progress tracking; early division of labor among team members (data, modeling, validation, documentation). |

**Timeline**

| Week | Tasks & Milestones |
| ----- | ----- |
| 1 | **Project Setup and Dataset Preparation** — Review the Ye et al. (2025) FLP database and key CO₂ activation literature. Define project scope, metrics, and deliverables. Extract and clean the CO₂-focused subset from FLPDB, removing incomplete or duplicate entries. Prepare a small validation subset for DFT benchmarking. |
| 2 | **Feature Engineering and Baseline Modeling** — Compute essential molecular descriptors (steric, electronic, and topological) using RDKit and xTB tools. Train baseline regression models (Random Forest, XGBoost) to predict CO₂ activation or binding energies. Evaluate performance using MAE, RMSE, and R² metrics, and visualize structure–property trends. |
| 3 | **Model Refinement and Candidate Screening** — Optimize model hyperparameters through cross-validation. Use model uncertainty to identify 5–10 promising new FLP candidates. Implement a simple active-learning iteration to prioritize compounds for DFT validation. |
| 4 | **DFT Validation and Finalization** — Perform DFT calculations on top candidates to confirm predicted CO₂ activation energies and compare results with ML outputs. Compile the curated dataset, trained models, and analysis. Prepare the final report and presentation, and release all code and data under an open-source license. |

**References:**

1. Ye, J., Wijethunga, C., & McEwen, M. (2025). The First Frustrated Lewis Pairs Database: Machine Learning and Cheminformatics-Aided Prediction of Small Molecule Activation. *The Journal of Physical Chemistry C*, *129*, 14346–14355. (from attached file)  
2. Khan, Md. N., van Ingen, Y., Boruah, T., McLauchlan, A., Wirth, T., & Melen, R. L. (2023). Advances in CO2 activation by frustrated Lewis pairs: from stoichiometric to catalytic reactions. *Chemical Science*, *14*(47), 13661–13681.  
3. Stephan, D. W. (2015). Frustrated Lewis pairs: from concept to catalysis. *Accounts of Chemical Research*, *48*(2), 306–316.  
4. Mömming, C. M., Otten, E., Kehr, G., Fröhlich, R., Grimme, S., Stephan, D. W., & Erker, G. (2009). Reversible metal-free hydrogen activation: a stoichiometric and catalytic example. *Angewandte Chemie International Edition*, *48*(36), 6643–6646.  
5. Pressler, J. P., & Gunanathan, C. (2017). Frustrated Lewis pair-catalyzed hydrogenation of carbonyl compounds, imines, and enamines. *ACS Catalysis*, *7*(11), 7433–7445.  
6. Butler, K. T., Davies, D. W., Cartwright, H., Isayev, O., & Walsh, A. (2018). Machine learning for molecular and materials science. *Nature*, *559*(7715), 547–555.  
7. Schleder, G. R., Padilha, A. C., Acosta, C. M., Costa, M., & Fazzio, A. (2019). From DFT to machine learning: recent advances in materials science. *Journal of Physics: Materials*, *2*(3), 032001\.

# **Ethics and Responsible Use:**

This project will be conducted with a strong commitment to ethical research, transparency, and open science. All data and results will be carefully curated, documented, and shared through an open-access repository to ensure reproducibility and foster collaboration. We will rely on well-established open-source tools for machine learning and cheminformatics, and release our code and datasets under a permissive license to encourage broad adoption. Throughout the work, we will remain attentive to the potential dual-use implications of AI in chemistry and ensure that all outcomes are directed toward advancing scientific understanding, environmental sustainability, and societal benefit.

## 

