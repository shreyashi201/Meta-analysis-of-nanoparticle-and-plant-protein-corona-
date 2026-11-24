# Meta-Analysis-and-Machine-Learning-of-Nanoparticle-Plant-Protein-Corona
Computational analysis, meta-analysis, and machine learning of nanoparticle–plant protein corona interactions, including proteomic data processing, classification, and identification of discriminatory proteins.

## Repository Structure

The repository contains the original R scripts used for machine learning-based analysis of nanoparticle–plant protein corona proteomic profiles, along with the corresponding datasets and analytical outputs.

### 1. LightGBM Analysis

The original R script (LightGBM_Original.R) implements a multiclass classification model using the Light Gradient Boosting Machine (LightGBM) algorithm.

The model was developed to classify proteomic profiles across eight experimental groups using repeated stratified cross-validation (10 repeats × 3 folds).

The analysis includes:

- Data preprocessing and median imputation.

- Multiclass model development and evaluation.

- Performance assessment using accuracy, macro-precision, macro-recall, and macro-F1 score.

- Confusion matrix generation.

- Cross-validation-based protein feature importance analysis using gain, cover, and split frequency.

- Identification and ranking of proteins contributing to classification.

The corresponding datasets, model performance metrics, cross-validation predictions, confusion matrices, and protein importance outputs are provided in the respective folders.

### 2. Random Forest Analysis

The original R script (RandomForest_Original.R) implements a multiclass classification approach using the Random Forest (RF) algorithm.

The model was developed using 1,000 decision trees and repeated stratified 3-fold cross-validation (10 repeats).

The analysis includes:

- Proteomic data preprocessing.

- Multiclass classification across eight experimental groups.

- Model performance evaluation using accuracy, precision, recall, and F1 score.

- Confusion matrix generation.

- Protein feature importance assessment using Mean Decrease in Gini.

- Ranking of proteins contributing to classification.

The corresponding datasets, model performance metrics, confusion matrices, and protein importance outputs are provided in the respective folders.

### 3. Support Vector Machine (SVM) Analysis

The original R script (SVM_Original.R) implements a multiclass Support Vector Machine (SVM) classification model using a linear kernel.

The model was evaluated using repeated stratified 3-fold cross-validation (10 repeats).

The analysis includes:

- Data preprocessing, median imputation, and feature standardization.

- Multiclass classification across eight experimental groups.

- Model performance evaluation using accuracy, macro-precision, macro-recall, and macro-F1 score.

- Confusion matrix generation.

- Evaluation of classification performance across repeated cross-validation runs.

The corresponding datasets, model performance metrics, and confusion matrix outputs are provided in the respective folders.

## Experimental Groups

The proteomic dataset comprises 24 samples representing eight experimental groups, with three biological replicates per group:

1. Fe3O4 Positive In vivo

2. Fe3O4 Positive In vitro

3. Fe3O4 Negative In vivo

4. Fe3O4 Negative In vitro

5. AuNP Green Hard Corona

6. AuNP Chemical Hard Corona

7. AuNP Green Soft Corona

8. AuNP Chemical Soft Corona

## Computational Environment

The analyses were performed in R using the respective machine learning and data visualization packages. The original R scripts are provided to document the computational workflows used in the study.

## Reproducibility and Data Availability

The original analysis scripts and their corresponding analytical outputs are provided to support transparency and reproducibility. The input datasets should be used in the format specified in the respective R scripts.

## Contributors

Research and computational analysis: Shreyashi Bansal, Department of Botany, University of Delhi

Research supervision: Prof. Renu Deswal, Department of Botany, University of Delhi



The manuscript describing the files in this repository is currently being written.
