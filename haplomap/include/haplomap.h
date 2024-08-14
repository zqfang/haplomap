//
// Created by Zhuoqing Fang on 7/11/20.
//
#pragma once
#ifndef HAPLOMAP_H
#define HAPLOMAP_H

#include <getopt.h>
#include <cstring>


using namespace std;
/// main entry point for subcommmands

/// convert vcf to NIEHS Format for eblocks input
int main_convert(int argc, char **argv);
/// convert ensembl-vep annotation for eblock input
int main_annot(int argc, char **argv);
/// finding haploblocks, see ehaploblocks.cpp
int main_eblocks(int argc, char **argv);
/// anova test, see quantTraitMap.cpp
int main_ghmap(int argc, char **argv);
/// genetic relation matrix, see pca.cpp
int main_eigen(int argc, char **argv);

#endif // HAPLOMAP_H
