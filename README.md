# mcp
Efficient exact algorithm for the Maximum Clique Problem in large sparse networks. Based on CUBIS graph decomposition and Bron–Kerbosch search.

# myMCP: A Fast Maximum Clique Algorithm for Large Sparse Networks

> This repository contains the official Python implementation of the algorithm proposed in our paper:

**Tianlong Fan**, Wenjun Jiang, Yi-Cheng Zhang, and Linyuan Lü,  
**"A fast maximum clique algorithm based on network decomposition for large sparse networks,"**  
*Chaos, Solitons and Fractals*, Vol. 194, 2025, 116239.  
[[Link to Paper](https://doi.org/10.1016/j.chaos.2025.116239)]

## 🔍 Description

This project provides a fast and scalable exact algorithm for solving the Maximum Clique Problem (MCP) in large sparse networks. It introduces a novel decomposition approach based on Complete-Upper-Bound-Induced Subgraphs (CUBIS), allowing significant pruning of the original graph.

The algorithm performs two main steps:
1. Constructs the first CUBIS from the top core layers to search for a candidate maximum clique.
2. Optionally constructs a second CUBIS if necessary to refine the result.

Experiments on 50 real-world networks (up to 20 million nodes) demonstrate near-linear time complexity.

## 📂 File Structure

- `myMCP.py`: Main implementation of the proposed MCP algorithm.
- `README.md`: Project overview and instructions.

## 🚀 Usage

1. Modify `NetworkAddress` in the script to point to your input edge list.
2. Set `LayerNum` to define how many top core layers to use in the first decomposition step.
3. Run:


