#!/bin/bash
deepwalk --format edgelist --representation-size 64 --walk-length 40 --number-walks 16 --workers 8 --input ./output/edgelist/Karate.edgelist --output ./data/interim/deepwalk_embeddings/2022-01-09/Karate.deepwalk.txt
deepwalk --format edgelist --representation-size 64 --walk-length 40 --number-walks 16 --workers 8 --input ./output/edgelist/Islands.edgelist --output ./data/interim/deepwalk_embeddings/2022-01-09/Islands.deepwalk.txt
deepwalk --format edgelist --representation-size 64 --walk-length 40 --number-walks 16 --workers 8 --input ./output/edgelist/Autophagy.edgelist --output ./data/interim/deepwalk_embeddings/2022-01-09/Autophagy.deepwalk.txt
deepwalk --format edgelist --representation-size 64 --walk-length 40 --number-walks 16 --workers 8 --input ./output/edgelist/DUB.edgelist --output ./data/interim/deepwalk_embeddings/2022-01-09/DUB.deepwalk.txt
deepwalk --format edgelist --representation-size 64 --walk-length 40 --number-walks 16 --workers 8 --input ./output/edgelist/Bioplex3.edgelist --output ./data/interim/deepwalk_embeddings/2022-01-09/Bioplex3.deepwalk.txt
deepwalk --format edgelist --representation-size 64 --walk-length 40 --number-walks 16 --workers 8 --input ./output/edgelist/HCT.edgelist --output ./data/interim/deepwalk_embeddings/2022-01-09/HCT.deepwalk.txt
deepwalk --format edgelist --representation-size 64 --walk-length 40 --number-walks 16 --workers 8 --input ./output/edgelist/Bioplex2.edgelist --output ./data/interim/deepwalk_embeddings/2022-01-09/Bioplex2.deepwalk.txt
deepwalk --format edgelist --representation-size 64 --walk-length 40 --number-walks 16 --workers 8 --input ./output/edgelist/Rep_Net.edgelist --output ./data/interim/deepwalk_embeddings/2022-01-09/Rep_Net.deepwalk.txt
