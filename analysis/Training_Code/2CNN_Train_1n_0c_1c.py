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

path = "./Train_Images"

# [:170000]
# [170000:195000]
if int(file_index) == 1:
    ggh_image_train = np.load(path + "/ggh_image_train_1n0c1c.npz")["arr_0"][:85000]
    ggh_jet_train = np.load(path + "/ggh_jet_train_1n0c1c.npz")["arr_0"][:85000]
    ggh_image_val = np.load(path + "/ggh_image_val_1n0c1c.npz")["arr_0"][:17500]
    ggh_jet_val = np.load(path + "/ggh_jet_val_1n0c1c.npz")["arr_0"][:17500]

    vbf_image_train = np.load(path + "/vbf_image_train_1n0c1c.npz")["arr_0"][:85000]
    vbf_jet_train = np.load(path + "/vbf_jet_train_1n0c1c.npz")["arr_0"][:85000]
    vbf_image_val = np.load(path + "/vbf_image_val_1n0c1c.npz")["arr_0"][:17500]
    vbf_jet_val = np.load(path + "/vbf_jet_val_1n0c1c.npz")["arr_0"][:17500]

    vh_image_train = np.load(path + "/vh_image_train_1n0c1c.npz")["arr_0"][:85000]
    vh_jet_train = np.load(path + "/vh_jet_train_1n0c1c.npz")["arr_0"][:85000]
    vh_image_val = np.load(path + "/vh_image_val_1n0c1c.npz")["arr_0"][:17500]
    vh_jet_val = np.load(path + "/vh_jet_val_1n0c1c.npz")["arr_0"][:17500]

    tth_image_train = np.load(path + "/tth_image_train_1n0c1c.npz")["arr_0"][:85000]
    tth_jet_train = np.load(path + "/tth_jet_train_1n0c1c.npz")["arr_0"][:85000]
    tth_image_val = np.load(path + "/tth_image_val_1n0c1c.npz")["arr_0"][:17500]
    tth_jet_val = np.load(path + "/tth_jet_val_1n0c1c.npz")["arr_0"][:17500]
    
if int(file_index) == 2:
    ggh_image_train = np.load(path + "/ggh_image_train_1n0c1c.npz")["arr_0"][85000:]
    ggh_jet_train = np.load(path + "/ggh_jet_train_1n0c1c.npz")["arr_0"][85000:]
    ggh_image_val = np.load(path + "/ggh_image_val_1n0c1c.npz")["arr_0"][17500:]
    ggh_jet_val = np.load(path + "/ggh_jet_val_1n0c1c.npz")["arr_0"][17500:]

    vbf_image_train = np.load(path + "/vbf_image_train_1n0c1c.npz")["arr_0"][85000:]
    vbf_jet_train = np.load(path + "/vbf_jet_train_1n0c1c.npz")["arr_0"][85000:]
    vbf_image_val = np.load(path + "/vbf_image_val_1n0c1c.npz")["arr_0"][17500:]
    vbf_jet_val = np.load(path + "/vbf_jet_val_1n0c1c.npz")["arr_0"][17500:]

    vh_image_train = np.load(path + "/vh_image_train_1n0c1c.npz")["arr_0"][85000:]
    vh_jet_train = np.load(path + "/vh_jet_train_1n0c1c.npz")["arr_0"][85000:]
    vh_image_val = np.load(path + "/vh_image_val_1n0c1c.npz")["arr_0"][17500:]
    vh_jet_val = np.load(path + "/vh_jet_val_1n0c1c.npz")["arr_0"][17500:]

    tth_image_train = np.load(path + "/tth_image_train_1n0c1c.npz")["arr_0"][85000:]
    tth_jet_train = np.load(path + "/tth_jet_train_1n0c1c.npz")["arr_0"][85000:]
    tth_image_val = np.load(path + "/tth_image_val_1n0c1c.npz")["arr_0"][17500:]
    tth_jet_val = np.load(path + "/tth_jet_val_1n0c1c.npz")["arr_0"][17500:]



print("ggH Training: Event Image {}, Jet Image {} ".format(ggh_image_train.shape,ggh_jet_train.shape))
print("ggH Val.: Event Image {}, Jet Image {} ".format(ggh_image_val.shape,ggh_jet_val.shape))
print("VBF Training: Event Image {}, Jet Image {} ".format(vbf_image_train.shape,vbf_jet_train.shape))
print("VBF Val.: Event Image {}, Jet Image {} ".format(vbf_image_val.shape,vbf_jet_val.shape))
print("VH Training: Event Image {}, Jet Image {} ".format(vh_image_train.shape,vh_jet_train.shape))
print("VH Val.: Event Image {}, Jet Image {} ".format(vh_image_val.shape,vh_jet_val.shape))
print("ttH Training: Event Image {}, Jet Image {} ".format(tth_image_train.shape,tth_jet_train.shape))
print("ttH Val.: Event Image {}, Jet Image {} ".format(tth_image_val.shape,tth_jet_val.shape))



x_train_eve = np.concatenate((ggh_image_train, vbf_image_train))
x_train_eve = np.concatenate((x_train_eve, vh_image_train))
x_train_eve = np.concatenate((x_train_eve, tth_image_train))
x_val_eve = np.concatenate((ggh_image_val, vbf_image_val))
x_val_eve = np.concatenate((x_val_eve, vh_image_val))
x_val_eve = np.concatenate((x_val_eve, tth_image_val))

x_train_jet = np.concatenate((ggh_jet_train, vbf_jet_train))
x_train_jet = np.concatenate((x_train_jet, vh_jet_train))
x_train_jet = np.concatenate((x_train_jet, tth_jet_train))
x_val_jet = np.concatenate((ggh_jet_val, vbf_jet_val))
x_val_jet = np.concatenate((x_val_jet, vh_jet_val))
x_val_jet = np.concatenate((x_val_jet, tth_jet_val))

# # sample_weight = np.concatenate((ggh_weight_list[HJ_train["eventindex"]], vbf_weight_list[VBF_train["eventindex"]]))
# # sample_weight = np.concatenate((sample_weight,vh_weight_list[VH_train["eventindex"]]))

y_train = np.concatenate((np.full(len(ggh_image_train), 0), np.full(len(vbf_image_train), 1)))
y_train = np.concatenate((y_train, np.full(len(vh_image_train), 2)))
y_train = np.concatenate((y_train, np.full(len(tth_image_train), 3)))

y_val = np.concatenate((np.full(len(ggh_image_val), 0), np.full(len(vbf_image_val), 1)))
y_val = np.concatenate((y_val, np.full(len(vh_image_val), 2)))
y_val = np.concatenate((y_val, np.full(len(tth_image_val), 3)))


y_train = to_categorical(y_train)
y_val = to_categorical(y_val)


print("Event: Training {}  Val {}".format(x_train_eve.shape,x_val_eve.shape))
print("Jet: Training {}  Val {}".format(x_train_jet.shape,x_val_jet.shape))
print("Target: Training {}  Val {}".format(y_train.shape,y_val.shape))




############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime consumption : {:.4f} min\033[0;m".format(totaltime/60.))


if int(file_index) == 1:

    ######################################################################################
    """
    Model
    """
    def return_pad_me(padding):
        def pad_me(x):
            #FRANK# x[:,:,:y,:] slice x off from y at the given axis.
            return(tf.concat((x,x[:,:,:padding,:]),2))
    #         return(tf.concat((2,x,x[:,:,:padding,:])))
        return(pad_me)


    input_shape = x_train_eve[0].shape
    print("image shape",input_shape)

    model_event = Sequential(name = 'Sequential_for_event')
    model_event.add(Lambda(return_pad_me(4),
                     input_shape=input_shape, name = 'event'))
    model_event.add(Conv2D(32, kernel_size=(5, 5), strides=(1, 1),
                     activation='relu',
                     data_format='channels_first', name = 'event_2D_1'))
    model_event.add(Lambda(return_pad_me(1),
                     input_shape=input_shape, name = 'event_padding_1'))
    model_event.add(MaxPooling2D(pool_size=(2, 2), strides=(2, 2),data_format='channels_first', name = 'event_MaxPooling_1'))
    model_event.add(Lambda(return_pad_me(4),input_shape=input_shape, name = 'event_padding_2'))
    model_event.add(Conv2D(64, (5, 5), activation='relu',data_format="channels_first", name = 'event_2D_2'))
    model_event.add(Lambda(return_pad_me(1),input_shape=input_shape, name = 'event_padding_3'))
    model_event.add(MaxPooling2D(pool_size=(2, 2),data_format="channels_first", name = 'event_MaxPooling_2'))
    model_event.add(Flatten(name = 'event_flatten'))
    model_event.add(Dense(300, activation='relu', name = 'event_dense_1'))

    model_event.add(Dropout(0.1))


    model_jet = Sequential(name = 'Sequential_for_jet')
    model_jet.add(Conv2D(32, kernel_size=(5, 5), strides=(1, 1),
                     activation='relu',
                    data_format='channels_first',input_shape=input_shape, name = 'jet'))
    model_jet.add(MaxPooling2D(pool_size=(2, 2), strides=(2, 2),data_format='channels_first', name = 'jet_MaxPooling_1'))
    model_jet.add(Conv2D(64, (5, 5), activation='relu',data_format='channels_first', name = 'jet_2D_1'))
    model_jet.add(MaxPooling2D(pool_size=(2, 2),data_format='channels_first', name = 'jet_MaxPooling_2'))
    model_jet.add(Flatten(name = 'jet_flatten'))
    model_jet.add(Dense(300, activation='relu', name = 'jet_dense_1'))

    model_jet.add(Dropout(0.1))


    mergedOut = Concatenate()([model_event.output,model_jet.output])
    # mergedOut = Dense(1, activation='sigmoid')(mergedOut)
    mergedOut = Dense(4, activation='softmax')(mergedOut)

    newModel = Model([model_event.input,model_jet.input], mergedOut,name = 'Combined')


    model_opt = keras.optimizers.Adadelta()

    newModel.compile(loss="categorical_crossentropy",#keras.losses.binary_crossentropy
                  optimizer=model_opt,
                  metrics=['accuracy'])
    newModel.summary()

if int(file_index) == 2:
    newModel = load_model("./Models/2CNN_1n0c1c_1.h5")



######################################################################################
# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
############################################################################################################################################################
check_list=[]
csv_logger = CSVLogger("./Models/training_log_1n0c1c_"+str(file_index)+".csv")
checkpoint = ModelCheckpoint(
                    filepath="./Models/checkmodel_1n0c1c_"+str(file_index)+".h5",
                    save_best_only=True,
                    verbose=1)
earlystopping = EarlyStopping(
                        monitor="loss",
                        min_delta=0.01,
                        patience=40,
                        verbose=1,
                        mode="auto",
                        baseline=None,
                        restore_best_weights=False,
                    )
check_list.append(checkpoint)
check_list.append(csv_logger)
check_list.append(earlystopping)
newModel.fit(  (x_train_eve, x_train_jet), y_train,
                validation_data = ((x_val_eve, x_val_jet), y_val),
                batch_size=512,
                epochs=150,
                shuffle=True,
                callbacks=check_list,
                verbose=1)

newModel.save("./Models/2CNN_1n0c1c_"+str(file_index)+".h5")

############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime Cost : {:.4f} min\033[0;m".format(totaltime/60.))

