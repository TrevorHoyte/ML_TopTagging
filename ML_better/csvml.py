import functools
import pandas as pd
import numpy as np
import tensorflow as tf
import os
import pathlib


train_file = pathlib.Path("/home/thoyte/dataset/ml/train_preprocessed.csv")
test_file = pathlib.Path("/home/thoyte/dataset/ml/test_preprocessed.csv")

#train_file = pathlib.Path("/home/thoyte/dataset/ml/vect3_shift.csv")
#test_file = pathlib.Path("/home/thoyte/dataset/ml/vect3_shift_test.csv")


df = pd.read_csv(train_file)
df_test = pd.read_csv(test_file)

#trim data keep only 100 const
#list_columns=[]

#for i in range(301,601):
# label=str(i)
# df.pop(label)
# df_test.pop(label)

print(df.head())
print(df.dtypes)

print("-------------------")

print(df_test.head())
print(df_test.dtypes)

target = df.pop('qcd=0, ttbar=1')
target_test = df_test.pop('qcd=0, ttbar=1')

dataset = tf.data.Dataset.from_tensor_slices((df.values, target.values))
dataset_test = tf.data.Dataset.from_tensor_slices((df_test.values, target_test.values))

for feat, targ in dataset.take(5):
  print ('Features: {}, Target: {}'.format(feat, targ))



train_dataset = dataset.shuffle(len(df)).batch(5000)
test_dataset = dataset_test.shuffle(len(df_test)).batch(100)

def get_compiled_model():
  model = tf.keras.Sequential([
    tf.keras.layers.Dense(300, activation='relu'),
    tf.keras.layers.Dense(102, activation='relu'),
    tf.keras.layers.Dense(12, activation='relu'),
    tf.keras.layers.Dense(6, activation='relu'),
    tf.keras.layers.Dense(1,activation='sigmoid')

    #tf.keras.layers.Dense(10, activation='relu'),
    #tf.keras.layers.Dense(10, activation='relu'),
    #tf.keras.layers.Dense(1)
  ])

  model.compile(optimizer='adam',
                loss=tf.keras.losses.BinaryCrossentropy(from_logits=True),
                metrics=['accuracy'])
  return model

model = get_compiled_model()
model.fit(train_dataset, epochs=40)

test_loss, test_accuracy = model.evaluate(test_dataset)

print('\n\nTest Loss {}, Test Accuracy {}'.format(test_loss, test_accuracy))

