#!/usr/bin/env python
# encoding: utf-8
from __future__ import absolute_import, division, print_function, unicode_literals
# Install TensorFlow
import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Dropout, Flatten , Convolution2D, MaxPooling2D , Lambda, Conv2D, Activation,Concatenate
from tensorflow.keras.optimizers import Adam , SGD , Adagrad
from tensorflow.keras.callbacks import ModelCheckpoint, LearningRateScheduler, EarlyStopping, CSVLogger, ReduceLROnPlateau
from tensorflow.keras.utils import to_categorical
from tensorflow.keras import regularizers , initializers
from tensorflow.keras.models import Model
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from tensorflow.keras.preprocessing.image import NumpyArrayIterator


gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
  # Restrict TensorFlow to only use the first GPU
    try:
        tf.config.experimental.set_visible_devices(gpus[0], 'GPU')
        tf.config.experimental.set_virtual_device_configuration(
        gpus[0],
        [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=12000)])
        logical_gpus = tf.config.experimental.list_logical_devices('GPU')
        print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPU")
    except RuntimeError as e:
    # Visible devices must be set before GPUs have been initialized
        print(e)
        
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.utils import shuffle
from sklearn.metrics import confusion_matrix
# from xgboost import XGBClassifier
import tensorflow.keras.backend as K
from sklearn import metrics

# Import local libraries
import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import importlib
import os
import sys

file_index = sys.argv[1]


print("Tensorflow Version is {}".format(tf.__version__))
print("Keras Version is {}".format(tf.keras.__version__))
from tensorflow.python.client import device_lib
# print(device_lib.list_local_devices())
tf.device('/device:XLA_GPU:0')

# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()

path = "./Test_Images"

if int(file_index) == 1:
    ggh_image_test = np.load(path + "/ggh_image_train_1n0c1c.npz")["arr_0"][:114500]
    ggh_jet_test = np.load(path + "/ggh_jet_train_1n0c1c.npz")["arr_0"][:114500]

    vbf_image_test = np.load(path + "/vbf_image_train_1n0c1c.npz")["arr_0"][:114500]
    vbf_jet_test = np.load(path + "/vbf_jet_train_1n0c1c.npz")["arr_0"][:114500]

    vh_image_test = np.load(path + "/vh_image_train_1n0c1c.npz")["arr_0"][:114500]
    vh_jet_test = np.load(path + "/vh_jet_train_1n0c1c.npz")["arr_0"][:114500]

    tth_image_test = np.load(path + "/tth_image_train_1n0c1c.npz")["arr_0"][:114500]
    tth_jet_test = np.load(path + "/tth_jet_train_1n0c1c.npz")["arr_0"][:114500]

if int(file_index) == 2:
    ggh_image_test = np.load(path + "/ggh_image_train_1n0c1c.npz")["arr_0"][114500:]
    ggh_jet_test = np.load(path + "/ggh_jet_train_1n0c1c.npz")["arr_0"][114500:]

    vbf_image_test = np.load(path + "/vbf_image_train_1n0c1c.npz")["arr_0"][114500:]
    vbf_jet_test = np.load(path + "/vbf_jet_train_1n0c1c.npz")["arr_0"][114500:]

    vh_image_test = np.load(path + "/vh_image_train_1n0c1c.npz")["arr_0"][114500:]
    vh_jet_test = np.load(path + "/vh_jet_train_1n0c1c.npz")["arr_0"][114500:]

    tth_image_test = np.load(path + "/tth_image_train_1n0c1c.npz")["arr_0"][114500:]
    tth_jet_test = np.load(path + "/tth_jet_train_1n0c1c.npz")["arr_0"][114500:]



print("ggH Test: Event Image {}, Jet Image {} ".format(ggh_image_test.shape,ggh_jet_test.shape))
print("VBF Test: Event Image {}, Jet Image {} ".format(vbf_image_test.shape,vbf_jet_test.shape))
print("VH Test: Event Image {}, Jet Image {} ".format(vh_image_test.shape,vh_jet_test.shape))
print("ttH Test: Event Image {}, Jet Image {} ".format(tth_image_test.shape,tth_jet_test.shape))



x_test_eve = np.concatenate((ggh_image_test, vbf_image_test))
x_test_eve = np.concatenate((x_test_eve, vh_image_test))
x_test_eve = np.concatenate((x_test_eve, tth_image_test))

x_test_jet = np.concatenate((ggh_jet_test, vbf_jet_test))
x_test_jet = np.concatenate((x_test_jet, vh_jet_test))
x_test_jet = np.concatenate((x_test_jet, tth_jet_test))


y_test = np.concatenate((np.full(len(ggh_image_test), 0), np.full(len(vbf_image_test), 1)))
y_test = np.concatenate((y_test, np.full(len(vh_image_test), 2)))
y_test = np.concatenate((y_test, np.full(len(tth_image_test), 3)))


y_test = to_categorical(y_test)


print("Event: Test {}".format(x_test_eve.shape))
print("Jet: Test {}".format(x_test_jet.shape))
print("Target: Test {}".format(y_test.shape))




############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime consumption : {:.4f} min\033[0;m".format(totaltime/60.))




######################################################################################
# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
############################################################################################################################################################


Model = load_model("./Models/2CNN_1n0c1c_2.h5")
prediction = Model.predict((x_test_eve, x_test_jet))
np.save("./Models/2CNN_1n_0c_1c_"+str(file_index),prediction)


############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime Cost : {:.4f} min\033[0;m".format(totaltime/60.))

