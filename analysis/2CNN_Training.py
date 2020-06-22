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


import csv_decoder
import save_and_load


print("Tensorflow Version is {}".format(tf.__version__))
print("Keras Version is {}".format(tf.keras.__version__))
from tensorflow.python.client import device_lib
print(device_lib.list_local_devices())
tf.device('/device:XLA_GPU:0')


######################################################################################
class DataGenerator:

    def __init__(self, dataframe, batch_size,):
  
        self.n_class = 4
        
        self.batch_size = batch_size
        self.total_len  = dataframe.shape[0] #training size
        self.dataframe = dataframe
        # self.rawX = dataX
        # self.rawY = dataY
        self.graddir_path_event = "/home/alan/ML_Analysis/higgs_classification/2CNN_Model/"
        self.graddir_path_jet = "/home/alan/ML_Analysis/higgs_classification/2CNN_Model/"
        self.on_epoch_end()

    def _build_pipeline(self, dataX1, dataX2, labelY):

        def loadnumpy_py(pathX1,pathX2):
            pathX1 = pathX1.numpy()
            pathX1 = os.path.join(self.graddir_path_event, pathX1.decode('UTF-8'))
            pathX2 = pathX2.numpy()
            pathX2 = os.path.join(self.graddir_path_jet, pathX2.decode('UTF-8'))

            grad1 = np.load(pathX1)["arr_0"]
            grad2 = np.load(pathX2)["arr_0"]
            return grad1, grad2
              

        def preprocess_fn(pathX1, pathX2, labelY):

            # processing fb
            labelY = tf.one_hot(labelY , self.n_class)
            
            # load image
            [dataX1, dataX2] = tf.py_function(loadnumpy_py, [pathX1, pathX2], [tf.float32,tf.float32])
            dataX1.set_shape([3, 40, 40])
            dataX2.set_shape([3, 40, 40])
            
            return (dataX1, dataX2), labelY

        dataset = tf.data.Dataset.from_tensor_slices( (dataX1, dataX2, labelY) )
        dataset = dataset.shuffle(self.batch_size * 8000)
        dataset = dataset.repeat()
        dataset = dataset.map(preprocess_fn, num_parallel_calls=tf.data.experimental.AUTOTUNE)
        dataset = dataset.batch(self.batch_size)
        dataset = dataset.prefetch(tf.data.experimental.AUTOTUNE)

        self.dataset = dataset
        self.tf_iter = iter(dataset)

    def  __len__(self):

        return self.total_len // self.batch_size

    def __getitem__(self, index):

        return next(self.tf_iter)

    def on_epoch_end(self):

        dataX1 = np.array(self.dataframe["EventImage"])
        dataX2 = np.array(self.dataframe["JetImage"])
        labelY = np.array(self.dataframe["Y"])

        self._build_pipeline(dataX1, dataX2 ,labelY)
        
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


input_shape = np.load("./2CNN_Model/JetTest_1/x_test_jet_1.npz")["arr_0"].shape
print(input_shape)
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
        
######################################################################################

dataframe_train = pd.read_csv("./2CNN_Model/Train_dict.csv")
dataframe_val = pd.read_csv("./2CNN_Model/Val_dict.csv")
print("# of Training Data: {} ".format(dataframe_train.shape[0]))
train_gen = DataGenerator(dataframe_train, 512)
valid_gen = DataGenerator(dataframe_val, 512)



######################################################################################
# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
############################################################################################################################################################




check_list=[]
csv_logger = CSVLogger('./2CNN_Model/training_log.csv')
checkpoint = ModelCheckpoint(
                    filepath='./2CNN_Model/checkmodel.h5',
                    save_best_only=True,
                    verbose=1)
check_list.append(checkpoint)
check_list.append(csv_logger)
newModel.fit(train_gen.dataset,
                    steps_per_epoch = train_gen.__len__(),
                    validation_data = valid_gen.dataset,
                    validation_steps= valid_gen.__len__(),
                    epochs=100,
                    callbacks=check_list,
                    verbose=1)

newModel.save("./2CNN_Model/model.h5")

############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime Cost : {:.4f} min\033[0;m".format(totaltime/60.))



