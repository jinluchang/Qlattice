#pragma once

#include <array>

QLAT_START_NAMESPACE

class Coordinate: public std::array<int, DIM>
{
public:
        inline Coordinate(){
                std::array<int, DIM>::array();
                std::array<int, DIM>::fill(0);
        }

        inline Coordinate(int first, int second, int third, int fourth){
                std::array<int, DIM>::array();
                int *p = std::array<int, DIM>::data();
                *p = first; *(p + 1) = second; *(p + 2) = third; *(p + 3) = fourth;
        }

        inline int product() const {
                int ret = 1;
                for(int i = 0; i < this->size(); i++){
                        ret *= std::array<int, DIM>::operator[](i);
                }
                return ret;
        }
};

QLAT_END_NAMESPACE


