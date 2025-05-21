import h5py
import pandas as pd
import numpy as np

filename = "bubble_data_cpp.h5"
dataset_name = "bubble_info_table"

with h5py.File(filename, 'r') as hf:
    structured_array = hf[dataset_name][:]

print("\nNumPy Structured Array record (first element):")
print(structured_array[0])
print("\nNumPy Structured Array dtypes:")
print(structured_array.dtype)

# Convert NumPy structured array to Pandas DataFrame
# Manual construction to handle multi-dimensional fields like 'location' robustly.
data_for_pandas = {}
for name in structured_array.dtype.names:
    column_data = structured_array[name]
    # structured_array.dtype.fields[name][0] is the dtype of the field itself.
    # .shape will be non-empty if the field elements are arrays (e.g., (3,) for location)
    if structured_array.dtype.fields[name][0].shape:
        # For 'location', column_data is (N,3). list(column_data) creates a
        # list of N arrays, each of shape (3,). This is what Pandas likes.
        data_for_pandas[name] = list(column_data)
    else:
        # For scalar fields, column_data is (N,), which is fine.
        data_for_pandas[name] = column_data

df_from_numpy = pd.DataFrame(data_for_pandas)
print(df_from_numpy.head(10))