import numpy as np
import pandas as pd
import pytest
from biological_staging.staging import coerce_stage_value, assign_stage_scalar
from biological_staging.stats import linear_contrast

def test_invalid_stage_is_not_rounded_down():
    assert np.isnan(coerce_stage_value(1.9))
    assert np.isnan(coerce_stage_value(np.inf))
    assert coerce_stage_value('2.0') == 2

def test_missing_biomarker_remains_missing():
    assert np.isnan(assign_stage_scalar(np.nan,10))
    assert np.isnan(assign_stage_scalar(10,np.nan))

def test_unknown_contrast_term_cannot_silently_be_zero():
    with pytest.raises(KeyError):
        linear_contrast(pd.Series([1.],index=['Intercept']),pd.DataFrame([[1.]],index=['Intercept'],columns=['Intercept']),{'missing_slope':1.})
