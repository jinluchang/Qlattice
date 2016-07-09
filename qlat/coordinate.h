#pragma once

#include <array>

QLAT_START_NAMESPACE

struct Coordinate: public std::array<int, DIM>
{
        Coordinate(){
                std::array<int, DIM>::fill(0);
        }

        Coordinate(int first, int second, int third, int fourth){
                int *p = data();
                p[0] = first;
                p[1] = second;
                p[2] = third;
                p[3] = fourth;
        }

        int product() const {
                int ret = 1;
                for(int i = 0; i < size(); i++){
                        ret *= operator[](i);
                }
                return ret;
        }
};

QLAT_END_NAMESPACE


