
from google.colab import files
uploaded = files.upload()

import pandas as pd
import gzip

# Load the .gz file
with gzip.open("variant_summary.txt.gz", 'rt') as f:
    df = pd.read_csv(f, sep='\t', low_memory=False)

# Show shape and first few rows
print("Data shape:", df.shape)
df.head()

print(df.columns)

# Optional: filter for human data only
df = df[df['Type'].notnull()]  # Usually already filtered to humans, but can check `Species` column if needed

# Keep only Pathogenic and Benign variants
df = df[df['ClinicalSignificance'].isin(['Pathogenic', 'Benign'])]

# Remove entries without a proper review status
df = df[df['ReviewStatus'].str.contains('criteria', na=False)]

print("Filtered data shape:", df.shape)

#Extract important biological features and the label
selected_columns = [
    'Type',
    'GeneSymbol',
    'ClinicalSignificance',
    'OriginSimple',
    'Chromosome',
    'Start',
    'ReferenceAllele',
    'AlternateAllele',
    'Name',
    'ReviewStatus'
]

df_filtered = df[selected_columns].copy()
df_filtered.head()

#Label encoding
df_filtered['label'] = df_filtered['ClinicalSignificance'].map({'Benign': 0, 'Pathogenic': 1})

# Check if variant is missense or nonsense (based on Name column)
df_filtered['is_missense'] = df_filtered['Name'].str.contains('missense', case=False, na=False).astype(int)
df_filtered['is_nonsense'] = df_filtered['Name'].str.contains('nonsense', case=False, na=False).astype(int)

# Add is_frameshift
df_filtered['is_frameshift'] = df_filtered['Name'].str.contains('fs', case=False, na=False).astype(int)

# Add is_splice_site
df_filtered['is_splice_site'] = df_filtered['Name'].str.contains(r'\+|-', regex=True, na=False).astype(int)

import pandas as pd
import re

import pandas as pd
import re

# Define mutation type based on keywords in the 'Name' column
def classify_mutation(name):
    if pd.isna(name):
        return "other"

    name_lower = name.lower()

    # Frameshift mutation
    if 'fs' in name_lower:
        return "frameshift"

    # Nonsense mutation (stop codon)
    if 'ter' in name_lower or '*' in name:
        return "nonsense"

    # Missense mutation: e.g., p.Arg330Met or p.Glu615Gly
    if re.search(r'p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}', name):
        return "missense"

    # Splice site mutation: e.g., c.25-2A>G or c.892+48G>A
    if re.search(r'c\.\d+[\+\-]\d+', name):
        return "splice_site"

    # Default fallback
    return "other"

# --- Load your data ---
# Example: df_filtered = pd.read_csv("your_file.csv")

# Apply the classification function
df_filtered['mutation_type'] = df_filtered['Name'].apply(classify_mutation)

# Create boolean columns for each mutation type
df_filtered['is_missense'] = df_filtered['mutation_type'] == 'missense'
df_filtered['is_nonsense'] = df_filtered['mutation_type'] == 'nonsense'
df_filtered['is_frameshift'] = df_filtered['mutation_type'] == 'frameshift'
df_filtered['is_splice_site'] = df_filtered['mutation_type'] == 'splice_site'

# (Optional) View mutation type distribution
print(df_filtered['mutation_type'].value_counts())

df_filtered.sample(5)

df_filtered.to_csv("clinvar_cleaned_variants.csv", index=False)
print("Saved cleaned dataset as 'clinvar_cleaned_variants.csv'")

from google.colab import files
# Download the CSV to your system
files.download("clinvar_cleaned_variants.csv")

"""REMOVE DUPLICATES"""

import pandas as pd

# Load the dataset
df = pd.read_csv("clinvar_cleaned_variants.csv")  # Replace with actual file name

# Display original number of rows
print("Original number of rows:", len(df))

# 1. Remove completely identical duplicate rows
df_unique = df.drop_duplicates()

# 2. (Recommended) Remove duplicates based on key columns
# These columns likely define unique variants
key_cols = ['Type', 'GeneSymbol', 'Chromosome', 'Start', 'Name']
df_unique = df_unique.drop_duplicates(subset=key_cols)

# Display updated number of rows
print("Number of rows after removing duplicates:", len(df_unique))

# Save the cleaned data (optional)
df_unique.to_csv("cleaned_variants.csv", index=False)

from google.colab import files
# Download the CSV to your system
files.download("cleaned_variants.csv")

from google.colab import files
uploaded = files.upload()

import pandas as pd
df_unique = pd.read_csv("cleaned_variants.csv")

"""Bar plot for mutation type"""

import matplotlib.pyplot as plt
import seaborn as sns

# Set seaborn style
sns.set(style="whitegrid")

# Countplot of mutation types
plt.figure(figsize=(8, 5))
sns.countplot(data=df_unique, x='mutation_type', order=df_unique['mutation_type'].value_counts().index)
plt.title('Distribution of Mutation Types')
plt.xlabel('Mutation Type')
plt.ylabel('Count')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

"""Bar Plot of Clinical Significance (label: 1=Pathogenic, 0=Benign)"""

# Custom color palette: 0 (Benign) = green, 1 (Pathogenic) = red
# Change the keys to strings '0' and '1' to match the likely data type in the 'label' column
palette = {'0': 'seagreen', '1': 'crimson'}

plt.figure(figsize=(6, 4))
sns.countplot(data=df_unique, x='label', palette=palette)
plt.title('Pathogenic vs Benign Variants')
plt.xlabel('Clinical Significance (0 = Benign, 1 = Pathogenic)')
plt.ylabel('Count')
plt.xticks([0, 1], ['Benign', 'Pathogenic'])
plt.tight_layout()
plt.show()

"""Heatmap: Mutation Type vs Clinical Significance"""

# Create a cross-tab
heatmap_data = pd.crosstab(df_unique['mutation_type'], df_unique['label'])

# Plot heatmap
plt.figure(figsize=(8, 5))
sns.heatmap(heatmap_data, annot=True, fmt='d', cmap='Blues')
plt.title('Mutation Type vs Clinical Significance')
plt.xlabel('Clinical Significance (0 = Benign, 1 = Pathogenic)')
plt.ylabel('Mutation Type')
plt.tight_layout()
plt.show()

"""Gene frequency plot (Top 10)"""

plt.figure(figsize=(10, 5))
top_genes = df_unique['GeneSymbol'].value_counts().nlargest(10)
sns.barplot(x=top_genes.index, y=top_genes.values, palette='viridis')
plt.title('Top 10 Genes with Most Variants')
plt.xlabel('Gene Symbol')
plt.ylabel('Number of Variants')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

from scipy.stats import chi2_contingency

# Create contingency table
contingency = pd.crosstab(df_unique['mutation_type'], df_unique['label'])

# Run chi-square test
chi2, p, dof, expected = chi2_contingency(contingency)

print("Chi-square Test Results:")
print(f"Chi2 Statistic: {chi2:.3f}")
print(f"p-value: {p:.4f}")
print(f"Degrees of Freedom: {dof}")
if p < 0.05:
    print("✅ Statistically significant association between mutation type and pathogenicity.")
else:
    print("❌ No significant association found.")

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from sklearn.preprocessing import LabelEncoder

# Load data
df = pd.read_csv("cleaned_variants.csv")  # Replace with actual file name

# Select features
features = ['Type', 'mutation_type', 'is_missense', 'is_nonsense', 'is_frameshift', 'is_splice_site', 'OriginSimple', 'Chromosome']
X = df[features]
y = df['label']

# Encode categorical features
for col in ['Type', 'mutation_type', 'OriginSimple', 'Chromosome']:
    le = LabelEncoder()
    X[col] = le.fit_transform(X[col].astype(str))

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

"""Random forest"""

# Train model
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Predict
y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)[:, 1]

# Evaluation
print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))
print("\nClassification Report:\n", classification_report(y_test, y_pred))
print("ROC-AUC Score:", roc_auc_score(y_test, y_prob))

"""Model training and evaluation"""

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from sklearn.model_selection import train_test_split

# Assuming your data is already preprocessed into X and y
# X = features, y = label column
#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Define models
models = {
    'Logistic Regression': LogisticRegression(max_iter=1000, random_state=42),
    'XGBoost': XGBClassifier(eval_metric='logloss', random_state=42)
}

# Train and evaluate each model
for name, model in models.items():
    print(f"\nModel: {name}")
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)[:, 1]

    print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))
    print("\nClassification Report:\n", classification_report(y_test, y_pred))
    print("ROC-AUC Score:", roc_auc_score(y_test, y_proba))

"""FEATURE IMPORTANCE - XGBoost"""

import matplotlib.pyplot as plt
from xgboost import plot_importance
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier

# Define models (re-defining for clarity in this block, but it's already in your notebook)
models = {
    'Logistic Regression': LogisticRegression(max_iter=1000, random_state=42),
    'XGBoost': XGBClassifier(eval_metric='logloss', random_state=42)
}

# Dictionary to store trained models
trained_models_dict = {}

# Train and evaluate each model (modified to store trained models)
for name, model in models.items():
    print(f"\nModel: {name}")
    model.fit(X_train, y_train)
    # Store the trained model in the dictionary
    trained_models_dict[name] = model

# Access the trained XGBoost model from the new dictionary
xgb_model = trained_models_dict['XGBoost']

# Plot top 15 features by importance (Gain)
plt.figure(figsize=(10, 6))
plot_importance(xgb_model,
                max_num_features=15,
                importance_type='gain',
                height=0.5,
                show_values=False)
plt.title('Top 15 Feature Importances - XGBoost (by Gain)')
plt.tight_layout()
plt.show()

"""Hypertuning"""

from xgboost import XGBClassifier
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
import numpy as np

# 1. Define parameter grid
param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [3, 5, 7, 10],
    'learning_rate': [0.01, 0.05, 0.1, 0.2],
    'subsample': [0.6, 0.8, 1.0],
    'colsample_bytree': [0.6, 0.8, 1.0],
    'gamma': [0, 0.1, 0.2, 0.3],
    'reg_alpha': [0, 0.1, 1],
    'reg_lambda': [1, 1.5, 2]
}

# 2. Initialize model
xgb_model = XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)

# 3. Randomized Search with Cross-Validation
random_search = RandomizedSearchCV(
    estimator=xgb_model,
    param_distributions=param_grid,
    n_iter=30,  # Try 30 random combinations
    scoring='roc_auc',
    cv=3,
    verbose=2,
    random_state=42,
    n_jobs=-1
)

# 4. Fit to training data
random_search.fit(X_train, y_train)

# 5. Best model
best_xgb = random_search.best_estimator_
print("Best Parameters:\n", random_search.best_params_)

# 6. Evaluate
y_pred = best_xgb.predict(X_test)
y_proba = best_xgb.predict_proba(X_test)[:, 1]

print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))
print("\nClassification Report:\n", classification_report(y_test, y_pred))
print("ROC-AUC Score:", roc_auc_score(y_test, y_proba))

"""SHAP ANALYSIS"""

import shap
import matplotlib.pyplot as plt
import xgboost as xgb

# Use a TreeExplainer for XGBoost
explainer = shap.TreeExplainer(best_xgb)  # Replace with your tuned XGBoost model

# Calculate SHAP values for test data
shap_values = explainer.shap_values(X_test)

# Plot summary (bar plot of global feature importance)
shap.summary_plot(shap_values, X_test)

import joblib
joblib.dump(best_xgb, "best_xgb_model.pkl")

import pandas as pd

importance_df = pd.DataFrame({
    'Feature': X_train.columns,
    'Importance': best_xgb.feature_importances_
}).sort_values(by="Importance", ascending=False)

importance_df.to_csv("xgb_feature_importance.csv", index=False)

from google.colab import files

# Download model file
files.download("best_xgb_model.pkl")

# Download feature importance CSV
files.download("xgb_feature_importance.csv")
