# built-in modules
import itertools
from typing import Any, Dict, List, Tuple

# third-party modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shap
from sklearn.model_selection import cross_validate, train_test_split
from sklearn.preprocessing import StandardScaler


def profile_xy_and_models(xy_to_try: Dict[str, Tuple[pd.DataFrame, pd.Series]],
                          models_to_try: Dict[str, Any],
                          scoring: str = 'r2', n_cv: int = 5,
                          x_preprocessing: List[str] = None,
                          y_preprocessing: List[str] = None,
                          verbose: bool = True) -> Tuple[pd.DataFrame, plt.Axes, plt.Axes]:
    """
    This mega-function enables the testing of multiple ML models on multiple Xy datasets for an
    initial understanding of model/feature performances. Multiple different normalizations for the
    X and y can be specified. This function uses a cross-validation approach, so no distinct
    validation set is created (just a training set). Assumes that the xy passed in have already
    been split for training (you can use create_train_and_lockbox_sets for this)

    Output includes plots summarizing the training and testing scores for each of the model/xy
    combinations, as well as a DataFrame containing those same results for further processing.

    :param Dict[str, Tuple[pd.DataFrame, pd.Series]] xy_to_try: a dictionary containing the names
        (keys) and Xy combinations (values, as tuples) to try modeling
    :param Dict[str, Any] models_to_try: the ML models (sklearn estimators) to try
    :param str scoring: the scoring metric to use; defaults to r2
    :param int n_cv: the number of cross validations to perform
    :param List[str] x_preprocessing: a list of pre-processing steps to apply to the feature matrix.
        Options include: 'standard' (standardizing/z-scoring)
    :param List[str] y_preprocessing: a list of pre-processing steps to apply to the targets.
        Options include: 'standard' (standardizing/z-scoring), 'log' (log+1 transform), and
        'standard_log' (apply log then standard transforms in succession)
    :param bool verbose: indicate if status updates should be printed at relevant moments
    """

    # prepare a DataFrame to store the results of the cross-validated comparisons
    profile_result_df = pd.DataFrame(
        columns=['model', 'xy', 'train_score', 'val_score']
    )

    # iterate through model/XY combinations and run a cross-validation in each case
    for model_name, model in models_to_try.items():
        for xy_name, (x_train_raw, y_train_raw) in xy_to_try.items():

            # normalize the X and y as requested by the options;
            x_list = [('X_raw', x_train_raw)]
            if x_preprocessing is not None:
                for x_norm_scheme in x_preprocessing:
                    if x_norm_scheme == 'standard':
                        x_train_standard = StandardScaler().fit_transform(x_train_raw)
                        x_list.append(('X_standard', x_train_standard))
                    # add more x_norm_scheme here
                    else:
                        raise ValueError(f'X normalization scheme {x_norm_scheme} not recognized')
            y_list = [('y_raw', y_train_raw)]
            if y_preprocessing is not None:
                for y_norm_scheme in y_preprocessing:
                    if y_norm_scheme == 'standard':
                        y_train_standard = StandardScaler().fit_transform(
                            y_train_raw.values.reshape(-1, 1)
                        )
                        y_list.append(('y_standard', y_train_standard))
                    elif y_norm_scheme == 'log':
                        y_train_log = np.log2(y_train_raw.values + 1)
                        y_list.append(('y_log', y_train_log))
                    elif y_norm_scheme == 'log10':
                        y_train_log = np.log10(y_train_raw.values + 1).reshape(-1,1)
                        y_list.append(('y_log10', y_train_log))
                    elif y_norm_scheme == 'standard_log':
                        y_train_standard_log = StandardScaler().fit_transform(
                            np.log2(y_train_raw.values + 1).reshape(-1, 1)
                        )
                        y_list.append(('y_standard_log', y_train_standard_log))
                    # add more y_norm_scheme here
                    else:
                        raise ValueError(f'y normalization scheme {y_norm_scheme} not recognized')

            xy_final = {}
            for (x_lab, x_mat), (y_lab, y_mat) in itertools.product(x_list, y_list):
                xy_lab_final = f'{xy_name}__{x_lab}__{y_lab}'
                xy_final[xy_lab_final] = (x_mat, y_mat)

            for xy_full_name, (X_to_use, y_to_use) in xy_final.items():
                _verbose_print(f'{model_name}: {xy_full_name}', verbose)

                cv_result = cross_validate(model, X_to_use, y=y_to_use, cv=n_cv, scoring=scoring,
                                           return_train_score=True, n_jobs=4)

                cv_result_df = pd.DataFrame(data={
                                                'model': model_name,
                                                'xy': xy_full_name,
                                                'train_score': cv_result['train_score'],
                                                'val_score': cv_result['test_score']
                                            })

                profile_result_df = profile_result_df.append(cv_result_df)

    # plot the training and testing results
    _, (ax_train, ax_val) = plt.subplots(2, 1, figsize=(15, 10))

    sns.boxplot(x='model', y='train_score', data=profile_result_df, hue='xy', dodge=True,
                fliersize=0, ax=ax_train)
    sns.swarmplot(x='model', y='train_score', data=profile_result_df, hue='xy', dodge=True,
                  color='black', ax=ax_train)
    handles_train, labels_train = ax_train.get_legend_handles_labels()
    ax_train.legend(handles_train[:int(len(handles_train)/2)],
                    labels_train[:int(len(labels_train)/2)], loc='center left',
                    bbox_to_anchor=(1.05, 0.5))
    ax_train.set_xlabel('')
    ax_train.set_ylabel(f'Training {scoring}', fontsize=13)
    ax_train.set_ylim(0, 1)

    sns.boxplot(x='model', y='val_score', data=profile_result_df, hue='xy', dodge=True,
                fliersize=0, ax=ax_val)
    sns.swarmplot(x='model', y='val_score', data=profile_result_df, hue='xy', dodge=True,
                  color='black', ax=ax_val)
    handles_val, labels_val = ax_val.get_legend_handles_labels()
    ax_val.legend(handles_val[:int(len(handles_val)/2)], labels_val[:int(len(labels_val)/2)],
                  loc='center left', bbox_to_anchor=(1.05, 0.5))
    ax_val.set_xlabel('')
    ax_val.set_ylabel(f'Validation {scoring}', fontsize=13)
    ax_val.set_ylim(0, 1)

    return profile_result_df, ax_train, ax_val


def feature_importance(x: pd.DataFrame, y: pd.Series, model, model_type: str = 'tree'):
    """
    Given a specific xy and a model, use Shapley values to visualize the feature importances. This
    function will automatically split the provided x and y into training and validation sets, a
    step required for using the Shapley values package.

    :param pd.DataFrame x: the X matrix
    :param pd.Series y: the targets to use for cross-validated training
    :param model: any sklearn model for which to determine feature importances
    :param str model_type: the type of model in use; options are 'tree', 'other'; this
        will influence the choice of explainer from the shap package that is used
    """

    x_train, x_val, y_train, y_val = train_test_split(x, y, test_size=0.2)
    fit_model = model.fit(x_train, y_train)

    if model_type == 'tree':
        explainer = shap.TreeExplainer(fit_model)
    else:
        explainer = shap.Explainer(fit_model)
    shap_values = explainer.shap_values(x_val)
    shap.summary_plot(shap_values, x_val)


def one_hot_feature_importance(x: pd.DataFrame, y: pd.Series, model, scoring: str = 'r2',
                               n_cv: int = 5, figsize=(15, 6)) -> plt.Axes:
    """
    Given an Xy/model set of choice, compute the feature importances. This function works for
    pure one-hot feature sets only. Currently only works with RandomForest models as well.
    Uses a 5-fold cross-validation to establish feature importance ranges.

    :param pd.DataFrame x: the X matrix
    :param pd.Series y: the targets to use for cross-validated training
    :param model: the (RandomForest) model to use
    :param str scoring: the scoring metric to use
    :param int n_cv: the number of cross-validations to perform
    :param Tuple[int, int] figsize: the size of the plot for matplotlib
    :return plt.Axes feat_imp_ax: the Axes containing the plot of one-hot feature importances
    """

    cv = cross_validate(model, x, y, cv=n_cv, scoring=scoring, return_estimator=True, n_jobs=4)

    estimators = cv['estimator']
    feat_imps = [est.feature_importances_ for est in estimators]

    feat_bps = []
    bases = []
    imps = []
    for feat, feat_imps in zip(x.columns, zip(*feat_imps)):
        bp, base = feat.split('_')
        feat_bps += [int(bp)] * len(feat_imps)
        bases += [base] * len(feat_imps)
        imps += list(feat_imps)

    feat_imp_cv_df = pd.DataFrame(data={'bp': feat_bps, 'imp': imps, 'base': bases})

    _, ax = plt.subplots(figsize=figsize)
    sns.lineplot(x='bp', y='imp', hue='base', hue_order=['G', 'C', 'T', 'A'],
                 data=feat_imp_cv_df, ax=ax)
    ax.set_xlabel('bp from TSS', fontsize=13)
    ax.set_ylabel('Feature Importance', fontsize=13)

    return ax


def create_train_and_lockbox_sets(
            xy_to_try: Dict[str, Tuple[pd.DataFrame, pd.Series]], random_state=None
        ) -> Tuple[Dict[str, Tuple[pd.DataFrame, pd.Series]],
                   Dict[str, Tuple[pd.DataFrame, pd.Series]]]:
    """
    Given a set of Xy combinations to try modeling, separate training/validation and lockbox sets

    :param Dict[str, Tuple[pd.DataFrame, pd.Series]] xy_to_try: the Xy combinations to split
    :param int random_state: set the random_state for sklearn's train_test_split (used under the
        hood) if a reproducible lockbox set is desired)
    :return Tuple[Dict[str, Tuple[pd.DataFrame, pd.Series]],
        Dict[str, Tuple[pd.DataFrame, pd.Series]] train_and_lockbox: the train and lockbox xy dicts
    """
    xy_train = {}
    xy_lockbox = {}

    for name, (X, y) in xy_to_try.items():
        x_train, x_lockbox, y_train, y_lockbox = train_test_split(X, y, test_size=0.1,
                                                                  random_state=random_state)
        xy_train[name] = (x_train, y_train)
        xy_lockbox[name] = (x_lockbox, y_lockbox)

    return xy_train, xy_lockbox


def _verbose_print(s: str, verbose: bool):
    """
    Print a string based on a Boolean toggle
    """
    if verbose:
        print(s)
