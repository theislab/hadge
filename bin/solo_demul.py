#!/usr/bin/env python
import scvi
import scanpy as sc
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Parser for SOLO - Doublet finding')
# arguments mainly for solo
parser.add_argument('--rna_matrix_dir', help='Input folder to RNA expression matrix in 10x format.')
# for solo.train()
parser.add_argument('--max_epochs', help='Number of epochs to train for', type=int, default=400)
parser.add_argument('--lr', help='Learning rate for optimization.', type=float, default=0.001)
parser.add_argument('--train_size', help='Size of training set in the range between 0 and 1.', type=float, default=0.9)
parser.add_argument('--validation_size',
                    help='Size of the test set. If None, defaults to 1 - train_size. If train_size + validation_size < 1, the remaining cells belong to a test set.',
                    type=float, default=None)
parser.add_argument('--batch_size', help='Minibatch size to use during training.', type=int, default=128)
parser.add_argument('--early_stopping', help='Adds callback for early stopping on validation_loss.', type=lambda x: (str(x).lower() in ['true']), default=True)
parser.add_argument('--early_stopping_patience', help='Number of times early stopping metric can not improve over early_stopping_min_delta.', type=int, default=30)
parser.add_argument('--early_stopping_min_delta', help='Threshold for counting an epoch towards patience train().', type=float, default=0.0)
# for solo.predict()
parser.add_argument('--soft', help='Return probabilities instead of class label', default=False, type=lambda x: (str(x).lower() in ['true']))
parser.add_argument('--include_simulated_doublets', help='Return probabilities for simulated doublets as well.',
                    type=lambda x: (str(x).lower() in ['true']), default=False)
parser.add_argument('--assignmentOutSolo', help='Name for the CSV file containing the output of Solo prediction',
                    default="solo_prediction", type=str)
parser.add_argument('--outputdir', help='Output directory')
args = parser.parse_args()

param_list = [['rna_matrix_dir', args.rna_matrix_dir], ['max_epochs', args.max_epochs], ['lr', args.lr], ['train_size', args.train_size], ['validation_size', args.validation_size], ['batch_size', args.batch_size], ['early_stopping', args.early_stopping], ['early_stopping_patience', args.early_stopping_patience], ['early_stopping_min_delta', args.early_stopping_min_delta], ['soft', args.soft], ['include_simulated_doublets', args.include_simulated_doublets]]
 
param_df = pd.DataFrame(param_list, columns=['Argument', 'Value'])

if __name__ == '__main__':
    adata = sc.read_10x_mtx(args.rna_matrix_dir)
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train(max_epochs=args.max_epochs, lr=args.lr, train_size=args.train_size,
               validation_size=args.validation_size, batch_size=args.batch_size,
               early_stopping=args.early_stopping,
               early_stopping_patience=args.early_stopping_patience,
               early_stopping_min_delta=args.early_stopping_min_delta)
    prediction = solo.predict(args.soft, include_simulated_doublets=args.include_simulated_doublets)
    prediction.to_csv(args.outputdir + "/" + args.assignmentOutSolo + "_res.csv")
    param_df.to_csv(args.outputdir + "/params.csv", index=False)



