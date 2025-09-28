import sys
import numpy as np
from sklearn.model_selection import StratifiedKFold
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten
from tensorflow.keras.layers import Conv1D, MaxPooling1D
from tensorflow.keras.utils import to_categorical

nome_train = sys.argv[1].split(".")[0]

def load_data(file):
    with open(file, "r") as f:
        records = f.readlines()[1:]  

    classes = list(set(line.strip().split(",")[-1] for line in records))

    X, Y = [], []
    for seq in records:
        elements = seq.strip().split(",")
        features = elements[1:-1]
        label = elements[-1]
        X.append(features)
        Y.append(classes.index(label))

    X = np.array(X, dtype=float)
    Y = np.array(Y, dtype=int)

    # normalize
    data_max = np.amax(X)
    X = X / data_max if data_max != 0 else X

    return X, Y, len(classes), X.shape[1]

def create_model(nb_classes, input_length):
    model = Sequential()
    model.add(Conv1D(filters=5, kernel_size=5, padding='valid',
                     input_shape=(input_length, 1)))
    model.add(Activation('relu'))
    model.add(MaxPooling1D(pool_size=2, padding='valid'))

    model.add(Conv1D(filters=10, kernel_size=5, padding='valid'))
    model.add(Activation('relu'))
    model.add(MaxPooling1D(pool_size=2, padding='valid'))

    model.add(Flatten())

    # MLP
    model.add(Dense(500, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(nb_classes, activation='softmax'))

    model.compile(optimizer='adam',
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])
    return model

def train_and_evaluate_model(model, datatr, labelstr, datate, labelste, nb_classes):
    datatr = datatr.reshape((datatr.shape[0], datatr.shape[1], 1))
    labelstr = to_categorical(labelstr, nb_classes)
    labelste_bin = to_categorical(labelste, nb_classes)

    model.fit(datatr, labelstr, epochs=100, batch_size=20, verbose=0)

    datate = datate.reshape((datate.shape[0], datate.shape[1], 1))

    preds = np.argmax(model.predict(datate, verbose=0), axis=1)

    return preds, labelste

if __name__ == "__main__":
    n_folds = 10
    X, Y, nb_classes, input_length = load_data(sys.argv[1])

    i = 1
    kfold = StratifiedKFold(n_splits=n_folds, shuffle=True)
    for train, test in kfold.split(X, Y):
        model = create_model(nb_classes, input_length)
        pred, Y_test = train_and_evaluate_model(model, X[train], Y[train], X[test], Y[test], nb_classes)

        np.save(f"./results/preds_{nome_train}_{i}_tanh.npy", pred)
        np.save(f"./results/test_{nome_train}_{i}_tanh.npy", Y_test)
        i += 1
