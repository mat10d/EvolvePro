#!/usr/bin/env python3
"""
EvolvePro 프레임워크 활용 개선 스크립트 (전체 데이터 활용)
소량 데이터(107개)를 모두 활용한 단백질 변이체 활성도 예측
- 능동 학습 제거, 전체 데이터로 학습
- Leave-One-Out 교차 검증으로 성능 평가
- 다중 모델 앙상블 적용
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path
import time
import warnings
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import Ridge, Lasso, ElasticNet
from sklearn.svm import SVR
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.base import clone
from Bio import SeqIO
import pickle
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

# EvolvePro 경로 추가
sys.path.append('/home/ubuntu/EvolvePro')

try:
    from evolvepro.src.data import process_variant
    EVOLVEPRO_AVAILABLE = True
except ImportError:
    print("⚠️  EvolvePro 모듈을 찾을 수 없습니다.")
    EVOLVEPRO_AVAILABLE = False


class FullDataTrainer:
    """전체 데이터를 활용한 학습 클래스"""

    def __init__(self, use_pca=True, n_components=40):
        self.use_pca = use_pca
        self.n_components = n_components
        self.models = {}
        self.best_model_name = None
        self.best_model = None
        self.ensemble_weights = {}

    def _build_pipeline(self, model, n_features):
        """스케일러 → (PCA) → 모델을 파이프라인으로 조합. LOO CV 내부에서도 각 fold마다 독립적으로 fit되어 데이터 누수를 방지한다."""
        steps = [('scaler', StandardScaler())]
        if self.use_pca and n_features > self.n_components:
            steps.append(('pca', PCA(n_components=self.n_components, random_state=42)))
        steps.append(('model', clone(model)))
        return Pipeline(steps)

    def create_models(self):
        """다양한 모델 생성"""
        print("  🤖 모델 생성 중...")

        models = {
            'Ridge_weak': Ridge(alpha=1.0, random_state=42),
            'Ridge_medium': Ridge(alpha=5.0, random_state=42),
            'Ridge_strong': Ridge(alpha=10.0, random_state=42),
            'Lasso': Lasso(alpha=0.5, random_state=42, max_iter=5000),
            'ElasticNet': ElasticNet(alpha=1.0, l1_ratio=0.5, random_state=42, max_iter=5000),
            'SVR': SVR(kernel='rbf', C=10, gamma='scale', epsilon=0.1),
            'RandomForest_small': RandomForestRegressor(
                n_estimators=50,
                max_depth=5,
                min_samples_split=3,
                min_samples_leaf=2,
                max_features='sqrt',
                random_state=42,
                n_jobs=-1
            ),
            'RandomForest_medium': RandomForestRegressor(
                n_estimators=100,
                max_depth=7,
                min_samples_split=2,
                min_samples_leaf=1,
                max_features='sqrt',
                random_state=42,
                n_jobs=-1
            ),
            'GradientBoosting': GradientBoostingRegressor(
                n_estimators=50,
                max_depth=3,
                learning_rate=0.1,
                min_samples_split=3,
                min_samples_leaf=2,
                random_state=42
            )
        }

        print(f"     총 {len(models)}개 모델 생성 완료")
        return models

    def evaluate_model_loo(self, model, X, y, model_name):
        """Leave-One-Out 교차 검증으로 모델 평가"""
        print(f"\n  📊 {model_name} 평가 중...")

        try:
            loo = LeaveOneOut()

            # LOO 교차 검증으로 예측값 얻기
            y_pred = cross_val_predict(model, X, y, cv=loo, n_jobs=-1)

            # 각종 지표 계산
            mse = mean_squared_error(y, y_pred)
            rmse = np.sqrt(mse)
            mae = mean_absolute_error(y, y_pred)
            r2 = r2_score(y, y_pred)

            # 상관계수 계산
            if len(set(y)) > 1 and len(set(y_pred)) > 1:
                spearman, _ = spearmanr(y, y_pred)
                pearson, _ = pearsonr(y, y_pred)
            else:
                spearman, pearson = 0.0, 0.0

            metrics = {
                'mse': mse,
                'rmse': rmse,
                'mae': mae,
                'r2': r2,
                'spearman': spearman,
                'pearson': pearson,
                'predictions': y_pred
            }

            print(f"     R² = {r2:.3f}, Spearman = {spearman:.3f}, MAE = {mae:.3f}")

            return metrics

        except Exception as e:
            print(f"     ❌ 평가 실패: {e}")
            return None

    def train_all_models(self, X, y):
        """모든 모델 학습 및 평가"""
        print("\n🚀 전체 데이터로 모델 학습 시작...")
        print(f"  📊 훈련 데이터: {X.shape}")

        n_features = X.shape[1]
        models = self.create_models()

        results = {}
        best_score = -np.inf

        print("\n" + "="*80)
        print("모델별 LOO 교차 검증 결과 (Pipeline: 각 fold 내부에서 독립적으로 전처리)")
        print("="*80)

        for name, model in models.items():
            # 각 모델마다 scaler+PCA+model을 하나의 파이프라인으로 묶어 LOO CV에 넘긴다.
            # cross_val_predict가 fold마다 파이프라인 전체를 fit하므로 데이터 누수가 없다.
            pipeline = self._build_pipeline(model, n_features)
            metrics = self.evaluate_model_loo(pipeline, X, y, name)

            if metrics is not None:
                results[name] = metrics
                score = metrics['spearman']
                if score > best_score:
                    best_score = score
                    self.best_model_name = name

        print("\n" + "="*80)
        print(f"🏆 최고 성능 모델: {self.best_model_name} (Spearman = {best_score:.3f})")
        print("="*80)

        # 최고 모델 파이프라인을 전체 데이터로 재학습
        best_pipeline = self._build_pipeline(models[self.best_model_name], n_features)
        best_pipeline.fit(X, y)
        self.best_model = best_pipeline

        # PCA 설명 분산 출력 (best 파이프라인 기준)
        if 'pca' in best_pipeline.named_steps:
            pca = best_pipeline.named_steps['pca']
            print(f"  🔍 PCA: {n_features} → {self.n_components} 차원")
            print(f"     설명된 분산 비율: {pca.explained_variance_ratio_.sum():.3f}")

        # 앙상블을 위한 모든 모델 파이프라인 학습
        print("\n  🔄 앙상블을 위한 모든 모델 학습...")
        self.models = {}

        for name, model in models.items():
            if results.get(name) is not None:
                pipeline = self._build_pipeline(model, n_features)
                pipeline.fit(X, y)
                self.models[name] = pipeline

                spearman = results[name]['spearman']
                self.ensemble_weights[name] = max(0, spearman)

        total_weight = sum(self.ensemble_weights.values())
        if total_weight > 0:
            self.ensemble_weights = {k: v/total_weight for k, v in self.ensemble_weights.items()}

        print(f"     {len(self.models)}개 모델 학습 완료")
        print(f"     앙상블 가중치: {self.ensemble_weights}")

        return results

    def predict_with_ensemble(self, X):
        """앙상블 예측 — 각 파이프라인이 내부적으로 전처리를 수행하므로 raw X를 그대로 전달"""
        predictions = []
        weights = []

        for name, pipeline in self.models.items():
            pred = pipeline.predict(X)
            weight = self.ensemble_weights.get(name, 0)

            if weight > 0:
                predictions.append(pred)
                weights.append(weight)

        if len(predictions) > 0:
            predictions = np.array(predictions)
            weights = np.array(weights)
            ensemble_pred = np.average(predictions, axis=0, weights=weights)
            pred_std = np.std(predictions, axis=0)
        else:
            ensemble_pred = self.best_model.predict(X)
            pred_std = np.zeros(len(ensemble_pred))

        return ensemble_pred, pred_std

    def predict_single_model(self, X):
        """최고 성능 단일 모델로 예측 — 파이프라인이 전처리를 내부에서 처리"""
        predictions = self.best_model.predict(X)

        # 불확실성 추정 (RandomForest인 경우 트리별 예측 분산 사용)
        final_estimator = self.best_model.named_steps['model']
        if hasattr(final_estimator, 'estimators_'):
            X_transformed = self.best_model[:-1].transform(X)
            tree_predictions = np.array([
                tree.predict(X_transformed) for tree in final_estimator.estimators_
            ])
            pred_std = np.std(tree_predictions, axis=0)
        else:
            pred_std = np.ones(len(predictions)) * np.std(predictions) * 0.1

        return predictions, pred_std


def load_data():
    """모든 데이터 로드"""
    print("🔍 데이터 로드 중...")

    # 파일 경로
    embedding_path = "/data/Evolvepro/output/plm/esm/DPE/DPE_esm2_t48_15B_UR50D.csv"
    exp_path = "/home/ubuntu/EvolvePro/data/exp/exp_data/DPE/rounds/DPE_conversion_data_1.xlsx"
    wt_fasta_path = "/home/ubuntu/EvolvePro/output/wt_fasta/DPE_WT.fasta"

    # 임베딩 데이터 로드
    print("  📊 임베딩 데이터 로드...")
    embeddings = pd.read_csv(embedding_path, index_col=0)

    print(f"     임베딩 데이터 크기: {embeddings.shape}")

    # WT 시퀀스 로드
    print("  🧬 WT 시퀀스 로드...")
    WT_sequence = str(SeqIO.read(wt_fasta_path, "fasta").seq)
    print(f"     WT 시퀀스 길이: {len(WT_sequence)}")

    # 실험 데이터 로드
    print("  📋 실험 데이터 로드...")
    exp_data = pd.read_excel(exp_path)
    exp_data.columns = ['variant', 'activity']

    # process_variant 적용
    exp_data["processed_variant"] = exp_data["variant"].apply(
        lambda x: process_variant(x, WT_sequence)
    )

    # 중복 제거
    exp_data = exp_data.drop_duplicates(subset=['processed_variant']).reset_index(drop=True)

    print(f"     실험 데이터 크기: {exp_data.shape}")
    print(f"     활성도 범위: {exp_data['activity'].min():.3f} ~ {exp_data['activity'].max():.3f}")

    return embeddings, exp_data, WT_sequence


def prepare_training_data(embeddings, exp_data):
    """학습 데이터 준비"""
    print("\n🔧 학습 데이터 준비...")

    # 임베딩과 실험 데이터 매칭
    embedding_variants = set(embeddings.index)
    exp_variants = set(exp_data['processed_variant'])

    # 매칭되는 변이체 찾기
    matched_variants = embedding_variants & exp_variants

    print(f"  📊 데이터 매칭 결과:")
    print(f"     임베딩 변이체 수: {len(embedding_variants)}")
    print(f"     실험 변이체 수: {len(exp_variants)}")
    print(f"     매칭되는 변이체 수: {len(matched_variants)}")

    if len(matched_variants) == 0:
        print("❌ 매칭되는 변이체가 없습니다!")
        return None, None

    # 매칭된 데이터만 선택
    matched_embeddings = embeddings.loc[list(matched_variants)]
    matched_exp_data = exp_data[exp_data['processed_variant'].isin(matched_variants)]

    # 인덱스 정렬하여 순서 맞추기
    matched_exp_data = matched_exp_data.set_index('processed_variant')
    matched_exp_data = matched_exp_data.loc[matched_embeddings.index]

    # 특성 데이터 (X)와 타겟 데이터 (y) 분리
    X = matched_embeddings.select_dtypes(include=[np.number])
    y = matched_exp_data['activity']

    print(f"  🎯 최종 학습 데이터:")
    print(f"     특성 데이터 (X): {X.shape}")
    print(f"     타겟 데이터 (y): {y.shape}")
    print(f"     타겟 통계: min={y.min():.3f}, max={y.max():.3f}, mean={y.mean():.3f}, std={y.std():.3f}")

    return X, y


def predict_all_variants(trainer, embeddings, exp_data, use_ensemble=True):
    """모든 변이체에 대한 예측"""
    print("\n🔮 전체 변이체 예측...")

    # 특성 데이터 준비
    X_all = embeddings.select_dtypes(include=[np.number])

    if use_ensemble:
        print("  🎯 앙상블 예측 사용")
        predictions, prediction_std = trainer.predict_with_ensemble(X_all)
    else:
        print("  🎯 단일 최고 모델 예측 사용")
        predictions, prediction_std = trainer.predict_single_model(X_all)

    # 결과 데이터프레임 생성
    results = pd.DataFrame({
        'variant': embeddings.index,
        'predicted_activity': predictions,
        'prediction_std': prediction_std,
        'confidence': 1 / (1 + prediction_std)
    })

    # 실험 데이터와 병합
    exp_data_indexed = exp_data.set_index('processed_variant')
    results = pd.merge(results, exp_data_indexed[['activity', 'variant']],
                      left_on='variant', right_index=True, how='left')

    # 컬럼 이름 정리
    if 'variant_y' in results.columns:
        results = results.rename(columns={'variant_x': 'processed_variant', 'variant_y': 'original_variant'})
    else:
        results = results.rename(columns={'variant': 'processed_variant'})

    if 'activity' in results.columns:
        results = results.rename(columns={'activity': 'actual_activity'})

    # 예측값 기준 내림차순 정렬
    results = results.sort_values('predicted_activity', ascending=False)

    print(f"  📊 예측 완료:")
    print(f"     전체 변이체 수: {len(results)}")
    print(f"     예측값 범위: {results['predicted_activity'].min():.3f} ~ {results['predicted_activity'].max():.3f}")

    return results


def plot_results(results, cv_results, trainer, output_dir, y=None):
    """결과 시각화"""
    print("\n📊 결과 시각화...")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # 1. 실험값 vs LOO 예측값 산점도
    # full-data 학습 모델의 in-sample 예측(results['predicted_activity'])이 아닌,
    # 각 샘플이 학습에서 제외된 fold에서 예측된 LOO 값을 사용해 정직한 성능 평가를 한다.
    best_name = trainer.best_model_name
    if y is not None and best_name in cv_results:
        y_true = np.array(y)
        y_loo_pred = cv_results[best_name]['predictions']

        plt.figure(figsize=(10, 8))
        plt.scatter(y_true, y_loo_pred, alpha=0.6, s=100, edgecolors='black')

        min_val = min(y_true.min(), y_loo_pred.min())
        max_val = max(y_true.max(), y_loo_pred.max())
        plt.plot([min_val, max_val], [min_val, max_val], 'r--', lw=2, label='Perfect prediction')

        spearman, _ = spearmanr(y_true, y_loo_pred)
        pearson, _ = pearsonr(y_true, y_loo_pred)
        r2 = r2_score(y_true, y_loo_pred)

        plt.xlabel('Actual Activity', fontsize=12)
        plt.ylabel('LOO Predicted Activity', fontsize=12)
        plt.title(
            f'LOO Cross-Validation: Predicted vs Actual ({best_name})\n'
            f'Spearman: {spearman:.3f}, Pearson: {pearson:.3f}, R²: {r2:.3f}',
            fontsize=14
        )
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path / 'prediction_scatter.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("  ✅ 산점도 저장 완료 (LOO 교차 검증 기반)")

    # 2. 모델 성능 비교
    if cv_results:
        model_names = list(cv_results.keys())
        spearmans = [cv_results[m]['spearman'] for m in model_names]
        r2s = [cv_results[m]['r2'] for m in model_names]
        maes = [cv_results[m]['mae'] for m in model_names]

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        # Spearman
        axes[0].barh(model_names, spearmans, color='steelblue')
        axes[0].set_xlabel('Spearman Correlation', fontsize=12)
        axes[0].set_title('Spearman Correlation by Model', fontsize=14)
        axes[0].axvline(x=0, color='black', linestyle='--', linewidth=0.8)
        axes[0].grid(True, alpha=0.3, axis='x')

        # R²
        axes[1].barh(model_names, r2s, color='coral')
        axes[1].set_xlabel('R² Score', fontsize=12)
        axes[1].set_title('R² Score by Model', fontsize=14)
        axes[1].axvline(x=0, color='black', linestyle='--', linewidth=0.8)
        axes[1].grid(True, alpha=0.3, axis='x')

        # MAE
        axes[2].barh(model_names, maes, color='lightgreen')
        axes[2].set_xlabel('MAE', fontsize=12)
        axes[2].set_title('MAE by Model', fontsize=14)
        axes[2].grid(True, alpha=0.3, axis='x')

        plt.tight_layout()
        plt.savefig(output_path / 'model_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("  ✅ 모델 비교 그래프 저장 완료")


def save_results(results, cv_results, trainer, output_dir):
    """결과 저장"""
    print(f"\n💾 결과 저장...")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # 1. 전체 예측 결과 저장
    results.to_csv(output_path / "full_data_predictions.csv", index=False)
    print(f"  ✅ 전체 예측 결과: {output_path / 'full_data_predictions.csv'}")

    # 2. 상위 변이체 저장
    top_variants = results.head(50)
    top_variants.to_csv(output_path / "top_50_variants.csv", index=False)
    print(f"  ✅ 상위 50 변이체: {output_path / 'top_50_variants.csv'}")

    # 3. 실험 검증 데이터 저장
    exp_results = results.dropna(subset=['actual_activity'])
    if len(exp_results) > 0:
        exp_results.to_csv(output_path / "experimental_validation.csv", index=False)
        print(f"  ✅ 실험 검증 데이터: {output_path / 'experimental_validation.csv'}")

    # 4. 모델 성능 비교 저장 ('predictions' 배열은 CSV에 직렬화할 수 없으므로 제외)
    if cv_results:
        cv_df = pd.DataFrame({
            name: {k: v for k, v in metrics.items() if k != 'predictions'}
            for name, metrics in cv_results.items()
        }).T
        cv_df.to_csv(output_path / "model_performance.csv")
        print(f"  ✅ 모델 성능 비교: {output_path / 'model_performance.csv'}")

    # 5. 최고 모델 파이프라인 저장 (StandardScaler + PCA + 모델이 하나로 통합됨)
    try:
        with open(output_path / "best_model.pkl", "wb") as f:
            pickle.dump(trainer.best_model, f)
        print(f"  ✅ 최고 모델 파이프라인: {output_path / 'best_model.pkl'}")
    except Exception as e:
        print(f"  ⚠️ 모델 저장 오류: {e}")

    # 6. 앙상블 가중치 저장
    if trainer.ensemble_weights:
        weights_df = pd.DataFrame(list(trainer.ensemble_weights.items()),
                                 columns=['model', 'weight'])
        weights_df = weights_df.sort_values('weight', ascending=False)
        weights_df.to_csv(output_path / "ensemble_weights.csv", index=False)
        print(f"  ✅ 앙상블 가중치: {output_path / 'ensemble_weights.csv'}")


def main():
    """메인 실행 함수"""
    print("=" * 80)
    print("EvolvePro 전체 데이터 활용 버전")
    print("소량 데이터(107개) 전체를 활용한 단백질 변이체 활성도 예측")
    print("=" * 80)

    start_time = time.time()

    try:
        # 1. 데이터 로드
        embeddings, exp_data, WT_sequence = load_data()

        # 2. 학습 데이터 준비
        X, y = prepare_training_data(embeddings, exp_data)

        if X is None:
            print("❌ 학습 데이터 준비 실패")
            return

        # 3. 전체 데이터로 학습
        trainer = FullDataTrainer(use_pca=True, n_components=40)
        cv_results = trainer.train_all_models(X, y)

        # 4. 전체 변이체 예측
        results = predict_all_variants(trainer, embeddings, exp_data, use_ensemble=True)

        # 5. 결과 저장 및 시각화
        output_dir = "/data/Evolvepro/output/DPE_evolvepro_full_activity_data_PCA_40"
        save_results(results, cv_results, trainer, output_dir)
        plot_results(results, cv_results, trainer, output_dir, y)

        # 6. 최종 요약
        print("\n" + "="*80)
        print("📊 최종 결과 요약")
        print("="*80)
        print(f"  🏆 최고 모델: {trainer.best_model_name}")
        print(f"  📊 총 변이체 수: {len(results)}")
        print(f"  🎯 실험 검증 변이체 수: {len(results.dropna(subset=['actual_activity']))}")

        # 최고 모델의 LOO 교차 검증 성능 (정직한 평가)
        best_name = trainer.best_model_name
        if best_name in cv_results:
            loo_metrics = cv_results[best_name]
            print(f"  📈 LOO 교차 검증 성능 ({best_name}):")
            print(f"     - Spearman: {loo_metrics['spearman']:.3f}")
            print(f"     - R²:       {loo_metrics['r2']:.3f}")
            print(f"     - MAE:      {loo_metrics['mae']:.3f}")

    except Exception as e:
        print(f"❌ 실행 중 오류 발생: {e}")
        import traceback
        traceback.print_exc()

    elapsed = time.time() - start_time
    print("=" * 80)
    print(f"✅ 전체 실행 완료! 소요 시간: {elapsed:.2f}초")
    print("=" * 80)


if __name__ == "__main__":
    main()
