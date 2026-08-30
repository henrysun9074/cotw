import os
import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from skopt import BayesSearchCV
from skopt.space import Real, Integer, Categorical
import matplotlib.pyplot as plt
import seaborn as sns
import joblib

print("loaded packages\n")

coral_df = pd.read_csv("/work/hs325/World_Corals/Cleaned data CSVs/qc_data.csv")
met_df = pd.read_csv("/work/hs325/World_Corals/Cleaned data CSVs/merged_met_plot_df.csv")

coral_df['scleractinia'] = np.where(coral_df['host_order'] == 'Scleractinia', 1, 0)
met_df['refined_origin'] = met_df['refined_origin'].str.replace('Host', 'Coral')

met_df_coral = met_df

print("loaded and cleaned data\n")

X = coral_df[coral_df.columns.intersection(met_df_coral['metabolite'])]
y = coral_df['scleractinia']
X = X.to_numpy()
y = y.to_numpy()

scaler = StandardScaler()
X = scaler.fit_transform(X)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123, stratify=y)
X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.2, random_state=123, stratify=y_train)

# hpt search spaces
search_spaces = {
    "RF": (
        RandomForestClassifier(n_jobs=-1, random_state=123),
        {
            "n_estimators": Integer(100, 2000),
            "max_depth": Integer(3, 50),
            "max_features": Categorical(["sqrt", "log2"]),
            "min_samples_split": Integer(2, 20),
            "min_samples_leaf": Integer(1, 10),
        },
    ),
    "GB": (
        XGBClassifier(tree_method="hist", eval_metric="logloss", use_label_encoder=False, n_jobs=-1, random_state=123),
        {
            "n_estimators": Integer(100, 2000),
            "max_depth": Integer(3, 15),
            "learning_rate": Real(1e-3, 0.3, prior="log-uniform"),
            "subsample": Real(0.5, 1.0),
            "colsample_bytree": Real(0.5, 1.0),
            "min_child_weight": Integer(1, 10),
            "gamma": Real(0, 5),
        },
    ),
}

# hpt using 100-iteration bayesian opt
results = {}
for name, (estimator, space) in search_spaces.items():
    print(f"starting hpt tuning for {name}")
    bayes = BayesSearchCV(
        estimator=estimator,
        search_spaces=space,
        n_iter=100,
        scoring="accuracy",
        cv=5,
        n_jobs=-1,
        refit=True,
        random_state=123,
        verbose=0,
    )
    bayes.fit(X_train, y_train)
    print(f"done with {name}, best accuracy = {bayes.best_score_:.4f}")
    print("best params = ", bayes.best_params_)

    best = bayes.best_estimator_
    val_acc = accuracy_score(y_val, best.predict(X_val))
    print(f"{name} validation accuracy: {val_acc:.4f}\n")

    results[name] = {
        "bayes": bayes,
        "best_estimator": best,
        "val_acc": val_acc,
        "best_params": bayes.best_params_,
    }

os.makedirs("models", exist_ok=True)

for name, info in results.items():
    best_model = info["best_estimator"]
    print(f"\nevaluating {name} on test set")

    y_pred = best_model.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    cm = confusion_matrix(y_test, y_pred)

    print(f"{name} accuracy {acc:.4f}")
    print(classification_report(y_test, y_pred))

    # Plot confusion matrix
    plt.figure(figsize=(5, 4))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", cbar=False)
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.tight_layout()
    plt.show()

    # Save model 
    model_path = f"models/{name}_best_model.joblib"
    joblib.dump(best_model, model_path)
    print(f"saved {name} model to /models")

