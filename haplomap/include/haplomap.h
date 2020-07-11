//
// Created by Zhuoqing Fang on 7/11/20.
//
#pragma once
#ifndef HBCGM_HAPLOMAP_H
#define HBCGM_HAPLOMAP_H

#include <getopt.h>
#include <cstring>

#ifndef __HAPLOMAP_VER__
#define __HAPLOMAP_VER__ "0.1.0"
#endif

#if defined(__clang__)
#define __COMPILER__ "clang++"
#elif defined(WIN32)
#define ___COMPILER__ "MSVC"
#else
#define __COMPILER__ "g++"
#endif


using namespace std;

int main_eblocks(int argc, char **argv);
int main_ghmap(int argc, char **argv);

#endif //HBCGM_HAPLOMAP_H
