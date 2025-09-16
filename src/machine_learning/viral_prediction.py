"""Machine learning models for viral interaction prediction."""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, roc_curve
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import networkx as nx
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform
import shap
import warnings
warnings.filterwarnings('ignore')
from typing import Dict, List, Tuple, Optional, Union

from ..utils import get_logger, log_processing_step, get_config_value


class ViralInteractionPredictor:
    """Machine learning models for predicting viral interactions."""
    
    def __init__(self):
        """Initialize the viral interaction predictor."""
        self.logger = get_logger('machine_learning')
        
        # Load configuration
        self.test_size = get_config_value('machine_learning.test_size', 0.2)
        self.random_state = get_config_value('machine_learning.random_state', 42)
        self.cv_folds = get_config_value('machine_learning.cv_folds', 5)
        
        # Model storage
        self.models = {}
        self.scalers = {}
        self.feature_importance = {}
        self.shap_explainers = {}
        
        self.logger.info("Viral interaction predictor initialized")
    
    def prepare_interaction_features(self, 
                                   abundance_matrix: pd.DataFrame,
                                   metadata: pd.DataFrame = None) -> pd.DataFrame:
        """Prepare features for viral interaction prediction.
        
        Args:
            abundance_matrix: Viral abundance matrix
            metadata: Sample metadata
            
        Returns:
            Feature matrix for machine learning
        """
        self.logger.info("Preparing interaction features")
        
        features = pd.DataFrame(index=abundance_matrix.index)
        
        with log_processing_step("Basic abundance features", self.logger):
            # Raw abundances (log-transformed)
            log_abundance = np.log10(abundance_matrix + 1)
            for virus in abundance_matrix.columns:
                features[f'{virus}_log_abundance'] = log_abundance[virus]
            
            # Presence/absence
            presence = (abundance_matrix > 0).astype(int)
            for virus in abundance_matrix.columns:
                features[f'{virus}_present'] = presence[virus]
        
        with log_processing_step("Diversity features", self.logger):
            # Sample-level diversity metrics
            features['viral_richness'] = (abundance_matrix > 0).sum(axis=1)
            features['viral_diversity_shannon'] = self._calculate_shannon_diversity(abundance_matrix)
            features['viral_evenness'] = self._calculate_evenness(abundance_matrix)
            features['total_viral_load'] = abundance_matrix.sum(axis=1)
        
        with log_processing_step("Co-occurrence features", self.logger):
            # Pairwise co-occurrence features
            viruses = abundance_matrix.columns.tolist()
            for i, virus_a in enumerate(viruses):
                for virus_b in viruses[i+1:]:
                    # Co-occurrence indicator
                    cooccur = ((abundance_matrix[virus_a] > 0) & (abundance_matrix[virus_b] > 0)).astype(int)
                    features[f'{virus_a}_{virus_b}_cooccur'] = cooccur
                    
                    # Abundance ratio (when both present)
                    ratio = np.where(
                        cooccur == 1,
                        np.log10((abundance_matrix[virus_a] + 1) / (abundance_matrix[virus_b] + 1)),
                        0
                    )
                    features[f'{virus_a}_{virus_b}_ratio'] = ratio
        
        with log_processing_step("Temporal and spatial features", self.logger):
            if metadata is not None:
                # Add metadata features
                for col in metadata.columns:
                    if col == 'sample_id':
                        continue
                    
                    if metadata[col].dtype == 'object':
                        # Categorical variables
                        le = LabelEncoder()
                        features[f'meta_{col}'] = le.fit_transform(metadata[col].fillna('unknown'))
                    else:
                        # Numerical variables
                        features[f'meta_{col}'] = metadata[col].fillna(metadata[col].median())
                
                # Temporal features if date available
                if 'collection_date' in metadata.columns:
                    dates = pd.to_datetime(metadata['collection_date'])
                    features['meta_day_of_year'] = dates.dt.dayofyear
                    features['meta_month'] = dates.dt.month
                    features['meta_season'] = ((dates.dt.month % 12) // 3).map({0: 0, 1: 1, 2: 2, 3: 3})  # 0=winter, 1=spring, 2=summer, 3=fall
        
        self.logger.info(f"Prepared {features.shape[1]} features for {features.shape[0]} samples")
        return features
    
    def _calculate_shannon_diversity(self, abundance_matrix: pd.DataFrame) -> pd.Series:
        """Calculate Shannon diversity index for each sample."""
        def shannon_sample(row):
            abundances = row[row > 0]
            if len(abundances) == 0:
                return 0
            proportions = abundances / abundances.sum()
            return -np.sum(proportions * np.log(proportions))
        
        return abundance_matrix.apply(shannon_sample, axis=1)
    
    def _calculate_evenness(self, abundance_matrix: pd.DataFrame) -> pd.Series:
        """Calculate Pielou's evenness index for each sample."""
        def evenness_sample(row):
            abundances = row[row > 0]
            if len(abundances) <= 1:
                return 0
            proportions = abundances / abundances.sum()
            shannon = -np.sum(proportions * np.log(proportions))
            max_shannon = np.log(len(abundances))
            return shannon / max_shannon if max_shannon > 0 else 0
        
        return abundance_matrix.apply(evenness_sample, axis=1)
    
    def train_exclusion_classifier(self, 
                                 features: pd.DataFrame,
                                 exclusion_pairs: List[Tuple[str, str]],
                                 abundance_matrix: pd.DataFrame) -> Dict:
        """Train classifier to predict viral exclusion patterns.
        
        Args:
            features: Feature matrix
            exclusion_pairs: List of virus pairs showing exclusion
            abundance_matrix: Original abundance matrix
            
        Returns:
            Dictionary with training results
        """
        self.logger.info("Training exclusion classifier")
        
        # Prepare target variable
        y = self._create_exclusion_labels(exclusion_pairs, abundance_matrix)
        
        if y.sum() < 10:  # Need sufficient positive examples
            self.logger.warning("Insufficient exclusion examples for training")
            return {'error': 'Insufficient exclusion examples'}
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            features, y, test_size=self.test_size, 
            random_state=self.random_state, stratify=y
        )
        
        # Scale features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        results = {
            'model_performance': {},
            'feature_importance': {},
            'predictions': {},
            'training_info': {
                'n_samples': len(features),
                'n_features': features.shape[1],
                'n_exclusion_samples': y.sum(),
                'class_balance': y.value_counts().to_dict()
            }
        }
        
        # Train multiple models
        models_to_train = {
            'random_forest': RandomForestClassifier(
                n_estimators=100, random_state=self.random_state, class_weight='balanced'
            ),
            'gradient_boosting': GradientBoostingClassifier(
                n_estimators=100, random_state=self.random_state
            ),
            'logistic_regression': LogisticRegression(
                random_state=self.random_state, class_weight='balanced', max_iter=1000
            ),
            'svm': SVC(
                random_state=self.random_state, class_weight='balanced', probability=True
            )
        }
        
        best_model = None
        best_score = 0
        
        for model_name, model in models_to_train.items():
            try:
                with log_processing_step(f"Training {model_name}", self.logger):
                    # Use scaled features for SVM and logistic regression
                    if model_name in ['svm', 'logistic_regression']:
                        model.fit(X_train_scaled, y_train)
                        y_pred = model.predict(X_test_scaled)
                        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
                        
                        # Cross-validation
                        cv_scores = cross_val_score(model, X_train_scaled, y_train, cv=self.cv_folds)
                    else:
                        model.fit(X_train, y_train)
                        y_pred = model.predict(X_test)
                        y_pred_proba = model.predict_proba(X_test)[:, 1]
                        
                        # Cross-validation
                        cv_scores = cross_val_score(model, X_train, y_train, cv=self.cv_folds)
                    
                    # Evaluate model
                    auc_score = roc_auc_score(y_test, y_pred_proba)
                    
                    results['model_performance'][model_name] = {
                        'auc_score': auc_score,
                        'cv_mean': cv_scores.mean(),
                        'cv_std': cv_scores.std(),
                        'classification_report': classification_report(y_test, y_pred, output_dict=True)
                    }
                    
                    # Feature importance
                    if hasattr(model, 'feature_importances_'):
                        importance = pd.Series(
                            model.feature_importances_, 
                            index=features.columns
                        ).sort_values(ascending=False)
                        results['feature_importance'][model_name] = importance.head(20).to_dict()
                    elif hasattr(model, 'coef_'):
                        importance = pd.Series(
                            np.abs(model.coef_[0]), 
                            index=features.columns
                        ).sort_values(ascending=False)
                        results['feature_importance'][model_name] = importance.head(20).to_dict()
                    
                    # Store best model
                    if auc_score > best_score:
                        best_score = auc_score
                        best_model = model_name
                        self.models['exclusion_classifier'] = model
                        if model_name in ['svm', 'logistic_regression']:
                            self.scalers['exclusion_classifier'] = scaler
                    
                    self.logger.info(f"{model_name} AUC: {auc_score:.3f}")
                    
            except Exception as e:
                self.logger.error(f"Error training {model_name}: {e}")
                results['model_performance'][model_name] = {'error': str(e)}
        
        results['best_model'] = best_model
        results['best_score'] = best_score
        
        # SHAP explanations for best model
        if best_model and best_model in ['random_forest', 'gradient_boosting']:
            try:
                with log_processing_step("Computing SHAP explanations", self.logger):
                    explainer = shap.TreeExplainer(self.models['exclusion_classifier'])
                    shap_values = explainer.shap_values(X_test[:100])  # Limit for performance
                    
                    if isinstance(shap_values, list):  # Binary classification
                        shap_values = shap_values[1]
                    
                    self.shap_explainers['exclusion_classifier'] = {
                        'explainer': explainer,
                        'shap_values': shap_values,
                        'test_data': X_test[:100]
                    }
                    
                    # Feature importance from SHAP
                    shap_importance = pd.Series(
                        np.abs(shap_values).mean(axis=0),
                        index=features.columns
                    ).sort_values(ascending=False)
                    
                    results['shap_importance'] = shap_importance.head(20).to_dict()
                    
            except Exception as e:
                self.logger.warning(f"SHAP explanation failed: {e}")
        
        return results
    
    def _create_exclusion_labels(self, 
                               exclusion_pairs: List[Tuple[str, str]], 
                               abundance_matrix: pd.DataFrame) -> pd.Series:
        """Create binary labels for exclusion prediction.
        
        Args:
            exclusion_pairs: List of virus pairs showing exclusion
            abundance_matrix: Abundance matrix
            
        Returns:
            Binary series indicating exclusion samples
        """
        labels = pd.Series(0, index=abundance_matrix.index)
        
        for virus_a, virus_b in exclusion_pairs:
            if virus_a in abundance_matrix.columns and virus_b in abundance_matrix.columns:
                # Samples where exclusion pattern is observed
                # (one virus present, the other absent)
                exclusion_mask = (
                    ((abundance_matrix[virus_a] > 0) & (abundance_matrix[virus_b] == 0)) |
                    ((abundance_matrix[virus_a] == 0) & (abundance_matrix[virus_b] > 0))
                )
                labels[exclusion_mask] = 1
        
        return labels
    
    def cluster_viral_profiles(self, 
                             abundance_matrix: pd.DataFrame,
                             method: str = 'kmeans',
                             n_clusters: int = None) -> Dict:
        """Cluster samples based on viral profiles.
        
        Args:
            abundance_matrix: Viral abundance matrix
            method: Clustering method ('kmeans', 'dbscan', 'hierarchical')
            n_clusters: Number of clusters (for kmeans and hierarchical)
            
        Returns:
            Dictionary with clustering results
        """
        self.logger.info(f"Clustering viral profiles using {method}")
        
        # Prepare data
        log_abundance = np.log10(abundance_matrix + 1)
        
        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(log_abundance)
        
        results = {
            'method': method,
            'cluster_labels': None,
            'cluster_centers': None,
            'silhouette_score': None,
            'cluster_stats': {},
            'dimensionality_reduction': {}
        }
        
        # Perform clustering
        if method == 'kmeans':
            if n_clusters is None:
                # Determine optimal number of clusters using elbow method
                n_clusters = self._find_optimal_clusters(X_scaled, max_k=min(10, len(X_scaled)//2))
            
            clusterer = KMeans(n_clusters=n_clusters, random_state=self.random_state)
            cluster_labels = clusterer.fit_predict(X_scaled)
            
            results['cluster_labels'] = cluster_labels
            results['cluster_centers'] = clusterer.cluster_centers_
            results['n_clusters'] = n_clusters
            
        elif method == 'dbscan':
            # Estimate eps using k-distance graph
            eps = self._estimate_dbscan_eps(X_scaled)
            
            clusterer = DBSCAN(eps=eps, min_samples=3)
            cluster_labels = clusterer.fit_predict(X_scaled)
            
            results['cluster_labels'] = cluster_labels
            results['eps'] = eps
            results['n_clusters'] = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
            
        elif method == 'hierarchical':
            if n_clusters is None:
                n_clusters = min(5, len(X_scaled)//10)
            
            # Compute linkage matrix
            linkage_matrix = linkage(X_scaled, method='ward')
            cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust') - 1
            
            results['cluster_labels'] = cluster_labels
            results['linkage_matrix'] = linkage_matrix
            results['n_clusters'] = n_clusters
        
        # Calculate silhouette score
        if len(set(cluster_labels)) > 1:
            from sklearn.metrics import silhouette_score
            results['silhouette_score'] = silhouette_score(X_scaled, cluster_labels)
        
        # Cluster statistics
        cluster_df = pd.DataFrame({
            'sample_id': abundance_matrix.index,
            'cluster': cluster_labels
        })
        
        for cluster_id in set(cluster_labels):
            if cluster_id == -1:  # Noise points in DBSCAN
                continue
            
            cluster_samples = abundance_matrix.index[cluster_labels == cluster_id]
            cluster_abundance = abundance_matrix.loc[cluster_samples]
            
            results['cluster_stats'][cluster_id] = {
                'n_samples': len(cluster_samples),
                'mean_viral_richness': (cluster_abundance > 0).sum(axis=1).mean(),
                'dominant_viruses': cluster_abundance.mean().sort_values(ascending=False).head(5).to_dict(),
                'prevalence_profile': (cluster_abundance > 0).mean().to_dict()
            }
        
        # Dimensionality reduction for visualization
        with log_processing_step("Dimensionality reduction", self.logger):
            # PCA
            pca = PCA(n_components=2, random_state=self.random_state)
            pca_coords = pca.fit_transform(X_scaled)
            
            results['dimensionality_reduction']['pca'] = {
                'coordinates': pca_coords,
                'explained_variance_ratio': pca.explained_variance_ratio_,
                'components': pca.components_
            }
            
            # t-SNE (if not too many samples)
            if len(X_scaled) <= 1000:
                tsne = TSNE(n_components=2, random_state=self.random_state, perplexity=min(30, len(X_scaled)//4))
                tsne_coords = tsne.fit_transform(X_scaled)
                
                results['dimensionality_reduction']['tsne'] = {
                    'coordinates': tsne_coords
                }
        
        return results
    
    def _find_optimal_clusters(self, X: np.ndarray, max_k: int = 10) -> int:
        """Find optimal number of clusters using elbow method."""
        inertias = []
        k_range = range(2, min(max_k + 1, len(X)))
        
        for k in k_range:
            kmeans = KMeans(n_clusters=k, random_state=self.random_state)
            kmeans.fit(X)
            inertias.append(kmeans.inertia_)
        
        # Find elbow point
        if len(inertias) < 2:
            return 2
        
        # Calculate rate of change
        rate_of_change = np.diff(inertias)
        elbow_point = np.argmax(rate_of_change) + 2  # +2 because we start from k=2
        
        return min(elbow_point, max_k)
    
    def _estimate_dbscan_eps(self, X: np.ndarray) -> float:
        """Estimate eps parameter for DBSCAN using k-distance graph."""
        from sklearn.neighbors import NearestNeighbors
        
        k = min(4, len(X) - 1)
        neighbors = NearestNeighbors(n_neighbors=k)
        neighbors.fit(X)
        distances, _ = neighbors.kneighbors(X)
        
        # Sort distances to k-th nearest neighbor
        k_distances = np.sort(distances[:, k-1])
        
        # Find knee point (simple heuristic)
        knee_point = len(k_distances) * 3 // 4
        eps = k_distances[knee_point]
        
        return eps
    
    def predict_viral_interactions(self, 
                                 features: pd.DataFrame,
                                 virus_pairs: List[Tuple[str, str]] = None) -> Dict:
        """Predict viral interactions for new samples.
        
        Args:
            features: Feature matrix for prediction
            virus_pairs: Specific virus pairs to predict (optional)
            
        Returns:
            Dictionary with predictions
        """
        if 'exclusion_classifier' not in self.models:
            raise ValueError("No trained model available. Train a model first.")
        
        model = self.models['exclusion_classifier']
        
        # Apply scaling if needed
        if 'exclusion_classifier' in self.scalers:
            X = self.scalers['exclusion_classifier'].transform(features)
        else:
            X = features.values
        
        # Make predictions
        predictions = model.predict(X)
        prediction_probabilities = model.predict_proba(X)[:, 1]
        
        results = {
            'predictions': predictions,
            'probabilities': prediction_probabilities,
            'sample_ids': features.index.tolist(),
            'high_risk_samples': features.index[prediction_probabilities > 0.7].tolist(),
            'prediction_summary': {
                'n_samples': len(predictions),
                'n_predicted_exclusions': predictions.sum(),
                'mean_probability': prediction_probabilities.mean(),
                'max_probability': prediction_probabilities.max()
            }
        }
        
        return results
    
    def build_viral_network_ml(self, 
                             abundance_matrix: pd.DataFrame,
                             method: str = 'correlation') -> nx.Graph:
        """Build viral interaction network using machine learning approaches.
        
        Args:
            abundance_matrix: Viral abundance matrix
            method: Network inference method ('correlation', 'mutual_info', 'graphical_lasso')
            
        Returns:
            NetworkX graph of viral interactions
        """
        self.logger.info(f"Building viral network using {method}")
        
        # Prepare data
        log_abundance = np.log10(abundance_matrix + 1)
        
        G = nx.Graph()
        viruses = abundance_matrix.columns.tolist()
        G.add_nodes_from(viruses)
        
        if method == 'correlation':
            # Pearson correlation network
            corr_matrix = log_abundance.corr()
            
            # Add edges for significant correlations
            threshold = 0.3
            for i, virus_a in enumerate(viruses):
                for j, virus_b in enumerate(viruses[i+1:], i+1):
                    corr = corr_matrix.loc[virus_a, virus_b]
                    if abs(corr) > threshold:
                        G.add_edge(virus_a, virus_b, 
                                 weight=abs(corr),
                                 correlation=corr,
                                 interaction_type='positive' if corr > 0 else 'negative')
        
        elif method == 'mutual_info':
            # Mutual information network
            from sklearn.feature_selection import mutual_info_regression
            
            threshold = 0.1
            for i, virus_a in enumerate(viruses):
                for virus_b in viruses[i+1:]:
                    mi = mutual_info_regression(
                        log_abundance[[virus_b]], 
                        log_abundance[virus_a],
                        random_state=self.random_state
                    )[0]
                    
                    if mi > threshold:
                        G.add_edge(virus_a, virus_b, 
                                 weight=mi,
                                 mutual_info=mi,
                                 interaction_type='association')
        
        elif method == 'graphical_lasso':
            # Graphical Lasso for sparse precision matrix
            from sklearn.covariance import GraphicalLassoCV
            
            try:
                gl = GraphicalLassoCV(cv=3, random_state=self.random_state)
                gl.fit(log_abundance)
                
                precision_matrix = gl.precision_
                
                # Add edges for non-zero precision matrix entries
                for i, virus_a in enumerate(viruses):
                    for j, virus_b in enumerate(viruses[i+1:], i+1):
                        precision = precision_matrix[i, j]
                        if abs(precision) > 1e-6:  # Threshold for numerical precision
                            G.add_edge(virus_a, virus_b,
                                     weight=abs(precision),
                                     precision=precision,
                                     interaction_type='conditional_dependence')
            
            except Exception as e:
                self.logger.warning(f"Graphical Lasso failed: {e}")
                return G
        
        # Add node attributes
        for virus in viruses:
            prevalence = (abundance_matrix[virus] > 0).mean()
            mean_abundance = abundance_matrix[virus].mean()
            
            G.nodes[virus].update({
                'prevalence': prevalence,
                'mean_abundance': mean_abundance,
                'degree': G.degree(virus)
            })
        
        self.logger.info(f"Built network with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
        
        return G
    
    def export_models(self, output_path: str) -> None:
        """Export trained models and results.
        
        Args:
            output_path: Output directory path
        """
        import pickle
        from pathlib import Path
        
        output_dir = Path(output_path)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save models
        for model_name, model in self.models.items():
            with open(output_dir / f'{model_name}.pkl', 'wb') as f:
                pickle.dump(model, f)
        
        # Save scalers
        for scaler_name, scaler in self.scalers.items():
            with open(output_dir / f'{scaler_name}_scaler.pkl', 'wb') as f:
                pickle.dump(scaler, f)
        
        # Save feature importance
        if self.feature_importance:
            importance_df = pd.DataFrame(self.feature_importance)
            importance_df.to_csv(output_dir / 'feature_importance.csv')
        
        self.logger.info(f"Models exported to {output_dir}")
    
    def load_models(self, input_path: str) -> None:
        """Load trained models.
        
        Args:
            input_path: Input directory path
        """
        import pickle
        from pathlib import Path
        
        input_dir = Path(input_path)
        
        # Load models
        for model_file in input_dir.glob('*.pkl'):
            if 'scaler' not in model_file.name:
                model_name = model_file.stem
                with open(model_file, 'rb') as f:
                    self.models[model_name] = pickle.load(f)
        
        # Load scalers
        for scaler_file in input_dir.glob('*_scaler.pkl'):
            scaler_name = scaler_file.stem.replace('_scaler', '')
            with open(scaler_file, 'rb') as f:
                self.scalers[scaler_name] = pickle.load(f)
        
        self.logger.info(f"Models loaded from {input_dir}")