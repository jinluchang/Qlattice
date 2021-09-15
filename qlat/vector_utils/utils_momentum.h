// utils_momentum.h
// Yi-Bo Yang
// Sept. 2014

#ifndef _UTILS_MOMENTUM_H_
#define _UTILS_MOMENTUM_H_

#pragma once

#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
///////#include "../kentucky/utils_lat_exchanger.h"

#ifndef _PI_
#define PI 3.1415926535898
#endif

#define _mom_order_ 100

namespace qlat{

//// layout.h
//// Ben Gamari
//// August 2009
struct dimensions {
        unsigned short nx, ny, nz, nt;

        dimensions() {}
        dimensions(unsigned short nx, unsigned short ny, unsigned short nz, unsigned short nt) :
                nx(nx), ny(ny), nz(nz), nt(nt) { }

        unsigned int compute_volume() const;
        unsigned short& operator[](unsigned int dim);
        unsigned short  operator[](unsigned int dim) const;
        dimensions operator*(dimensions& b) const;
};


struct momentum
{
    int iserial;
    int p2,p4;

    inline void ind_2_mom(int *momentum)
    {
      momentum[0]=sign_mom((iserial/(_mom_order_*_mom_order_))%_mom_order_);
      momentum[1]=sign_mom((iserial/(_mom_order_))%_mom_order_);
      momentum[2]=sign_mom(iserial%_mom_order_);
    }

    momentum(int iser=505050):iserial(iser)
    {
       int mom[3];
       ind_2_mom(mom);
       set(mom[0],mom[1],mom[2]);
    };

    momentum(int mom1,int mom2,int mom3)
    {
        int i11=mom1+_mom_order_/2;
        int i12=mom2+_mom_order_/2;
        int i13=mom3+_mom_order_/2;
        iserial=(i11*_mom_order_+i12)*_mom_order_+i13;
        set(mom1,mom2,mom3);
    }
    
    int operator[](int im)
    {
        if(im==0) return sign_mom((iserial/(_mom_order_*_mom_order_))%_mom_order_);
        if(im==1) return sign_mom((iserial/(_mom_order_))%_mom_order_);
        if(im==2) return sign_mom(iserial%_mom_order_);
        return 0;
    }

    momentum flip(){
    int mom[3];ind_2_mom(mom);
    return momentum(-mom[0],-mom[1],-mom[2]);
    }
    momentum plus(momentum b){
    int mom[6];
    ind_2_mom(mom);b.ind_2_mom(mom+3);
    for(int i=0;i<3;i++)mom[i]+=mom[i+3];
    return momentum(mom[0],mom[1],mom[2]);
    }
    momentum sub(momentum b){
    int mom[6];
    ind_2_mom(mom);b.ind_2_mom(mom+3);
    for(int i=0;i<3;i++)mom[i]-=mom[i+3];
    return momentum(mom[0],mom[1],mom[2]);
    }
    
    operator int(){ return this->iserial; }
    
    static int comparer(momentum &a,momentum &b)
    {
        if(a.p2>b.p2) return 1;
        else if(a.p2<b.p2) return -1;
        if(a.p4<b.p4) return 1;
        else if(a.p4>b.p4) return -1;
        return 0;
    }

private:

    inline void set(int mom1,int mom2,int mom3)
    {
       p2=mom1*mom1+mom2*mom2+mom3*mom3;
       p4=mom1*mom1*mom1*mom1
         +mom2*mom2*mom2*mom2
         +mom3*mom3*mom3*mom3;
    }

    inline int sign_mom(int ind)
    {
        if(ind>=20)return ind-_mom_order_/2;
        int mom=ind%(_mom_order_/10),sign=(ind/(_mom_order_/10)==0)?1:-1;
        return sign*mom;
    }

};

template<class T>
void generic_sort(T **list, int left,int right,int (*comparer)(T &_a,T &_b))
{
        int p = (left + right) / 2;
        ///T *tmp = list[p],*tmp2;
        T *tmp = list[p];
        int i = left,j = right;
        while (i < j) {
                while (i < p && comparer(*list[i],*tmp)<=0)
                        ++i;
                if (i < p) {
                        list[p] = list[i];
                        p = i;
                }
                while (j > p && comparer(*tmp,*list[j])<=0)
                        --j;
                if (j > p) {
                        list[p] = list[j];
                        p = j;
                }
        }
        list[p]=tmp;
        if(left<i)generic_sort(list,left, i - 1,comparer);
        if(right>i)generic_sort(list,i + 1, right,comparer);
}

template<class T>
void generic_sort(std::vector<T*> &list,int (*comparer)(T &_a,T &_b))
{
   generic_sort(list.data(),0,list.size()-1,comparer);
}

class mom_set
{
     int clear_flag;
     bool enhance_flag;

public:
    int nm;
    std::vector<int > mom_ind0,mode_count,mode_off;
    int * mom_ind;
    
    std::vector<std::vector<momentum> > mode;
    int n_mom,n_eff;

    int ind_2_mode(int ind){
    momentum indm(ind);
    for(int jq2=0;jq2<n_mom;jq2++)
    if(mode[jq2].size()>0)
    if(indm.p4==mode[jq2][0].p4 ||
       (enhance_flag==false && indm.p2==mode[jq2][0].p2) )
       return jq2;
    return 0;
    }
    int ind_2_serial(int ind){
    int imode=ind_2_mode(ind);
    for(unsigned int i=0;i<mode[imode].size();i++)
     if(mode[imode][i].iserial==ind) return mode_off[imode]+i;
    return 0;
    }

    inline void ind_2_mom(int ind,int *mom)
    {   momentum(ind).ind_2_mom(mom); }

    inline int mom_2_ind(int mom1,int mom2,int mom3)
    {   return momentum(mom1,mom2,mom3).iserial; }

    mom_set(int nm0=9):nm(nm0)
    {  
       n_mom=0;n_eff=0;clear_flag=0;enhance_flag=false;
    }
    
    void set(bool enhance=false){
    if(enhance==false)
    {
       mode.resize(3*nm*nm+1);
       for(int i=nm;i>=-nm;i--)
       for(int j=nm;j>=-nm;j--)
       for(int k=nm;k>=-nm;k--)
       {
          int ind=i*i+j*j+k*k;
          mode[ind].push_back(momentum(i,j,k));
       }
    }
    else
    {
       std::vector<momentum> lst;
       std::vector<momentum*> plst;
       for(int i=nm;i>=-nm;i--)
       for(int j=nm;j>=-nm;j--)
       for(int k=nm;k>=-nm;k--)
          lst.push_back(momentum(i,j,k));
       for(unsigned int i=0;i<lst.size();i++)
          plst.push_back(lst.data()+i);
       generic_sort(plst.data(),0,lst.size()-1,lst[0].comparer);
       int last_off=0;int present_size=0;
       for(unsigned int i=0;i<lst.size();i++)
       if(i==lst.size()-1 || lst[0].comparer(*plst[i+1],*plst[i])==1)
       {
           present_size++;
           mode.resize(present_size);
           for(unsigned int j=last_off;j<=i;j++)
              mode[present_size-1].push_back(*plst[j]);
          last_off=i+1;
       }
       enhance_flag=true;
    }
    
     int off=0;
     for(unsigned int i=0;i<mode.size();i++)
     {
          for(unsigned int j=0;j<mode[i].size();j++)mom_ind0.push_back(mode[i][j].iserial);
          mode_count.push_back(mode[i].size());
          mode_off.push_back(off);
          off+=mode[i].size();
          n_mom++;
         if(mode[i].size()>0)n_eff++;
     }
     mom_ind=mom_ind0.data();
    }

    int mom_flip(int ind1)
    { return momentum(ind1).flip().iserial; }

    int mom_add(int ind1,int ind2)
    { return momentum(ind1).plus(momentum(ind2)).iserial; }

    int mom_sub(int ind1,int ind2)
    { return momentum(ind1).sub(momentum(ind2)); }
    
    void print()
    {
       for(unsigned int i=0;i<mode.size();i++)
       if(mode[i].size()>0)
       {
          printf("%4d%6d%6d",i,mode[i][0].p2,mode[i][0].p4);
          for(unsigned int j=0;j<mode[i].size();j++)
             printf("%10d",mode[i][j].iserial);
          printf("\n");
       }
    }

};

struct mom_mode : momentum
{
private:    
    void set(int ix,int iy,int iz,double mass0,dimensions& desc)
    {
       mom[0]=2*sin(ix*PI/desc.nx);
       mom[1]=2*sin(iy*PI/desc.ny);
       mom[2]=2*sin(iz*PI/desc.nz);
       set_mass(mass0);
    }
public:
    double mom[4],mass;
    
    mom_mode(const mom_mode &a)
    {for(int i=0;i<4;i++)mom[i]=a.mom[i];mass=a.mass;iserial=a.iserial;}
    mom_mode(int ix,double mass0,dimensions& desc):momentum(ix)
    {
        int imom[3];ind_2_mom(imom);
        set(imom[0],imom[1],imom[2],mass0,desc);
    }
    
    mom_mode(int ix,int iy,int iz,double mass0,dimensions& desc):momentum(ix,iy,iz)
    {
         mom_set m_set(0);iserial=m_set.mom_2_ind(ix,iy,iz);
         set(ix,iy,iz,mass0,desc);
    }
    
    void set_mass(double mass0)
    {
       mass=mass0;
       mom[3]=sqrt(mass*mass+pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2));
    }
    
    void print(bool flag=false)
    {printf("(");for(int i=0;i<4;i++)printf("%10.5f",mom[i]);printf(")");
     if(flag)printf("\n");
    }
    
};

struct q2_mode
{
    mom_mode src,sink;
    double q2;
    q2_mode(mom_mode &a,mom_mode &b):src(a),sink(b)
    {
       q2=get_q2();
    }
    q2_mode(q2_mode &res):src(res.src),sink(res.sink),q2(res.q2){}
    
    double p2i(){return pow(src.mom[0],2)+pow(src.mom[1],2)+pow(src.mom[2],2);}
    double p2f(){return pow(sink.mom[0],2)+pow(sink.mom[1],2)+pow(sink.mom[2],2);}
    
    double get_q2(double mass1=-1.0,double mass2=-1.0){
    if(mass2<0&&mass1>=0)return get_q2(mass1,mass1);
    if(mass1>=0)
    {
       src.set_mass(mass1);sink.set_mass(mass2);
    }
    double  q20=0.0;
    for(int i=0;i<3;i++)
       q20+=pow(src.mom[i]-sink.mom[i],2);
    q20-=pow(src.mom[3]-sink.mom[3],2);
    return q20;
    }
    
    void print(bool flag=false)
    {
      printf("%13.6f: ",q2);
      src.print();printf(", ");sink.print();
      if(flag)printf("\n");
    }
};



struct q2_list
{
  std::vector<q2_mode *> list;
  int Nq2,size;
  std::vector<int> offset,Npmode;
  dimensions &desc;
  
  void set(std::vector<mom_mode> &src,std::vector<mom_mode> &sink,int cut_off=100,int iseed=0){
      qlat::RngState rs(abs(iseed)+1);
      ////qlat::randGen ran(rs);
      ///double_prn ran(abs(iseed)+1);
      mom_set momX;momX.set(true);
      offset.clear();Npmode.clear();
      for(unsigned int i=0;i<list.size();i++) delete list[i];
      list.clear();
      for(unsigned int i=0;i<src.size();i++)
      for(unsigned int j=0;j<sink.size();j++)
      if(momX.ind_2_mode(src[i].sub(sink[j]))<cut_off)
         list.push_back(new q2_mode(src[i],sink[j]));
      size=list.size();
      generic_sort(list.data(),0,size-1,this->comparer);
      
      int last_off=0;
      std::vector<q2_mode *> list_new;

      std::vector<int> sink_offset;
      if(iseed!=0)
      for(unsigned int i=0;i<momX.mode.size();i++)
      {
             //int ioff=(int)(ran.rand()*momX.mode[i].size());
             int ioff=(int)(qlat::u_rand_gen(rs)*momX.mode[i].size());
             sink_offset.push_back(ioff%momX.mode[i].size());
      }
      
      for(int i=0;i<size;i++)
      if(i==size-1||list[i+1]->q2-list[i]->q2>1e-7)
      {
         int size_t=i-last_off+1;
         if(iseed>0)
         {
             momentum iq=list[i]->src.sub(list[i]->sink);
             int iq2=momX.ind_2_mode(iq);
//             int ioff=(int)(ran.rand()*momX.mode[iq2].size());
//             momentum iqr=momX.mode[iq2][ioff%momX.mode[iq2].size()];
             momentum iqr=momX.mode[iq2][sink_offset[iq2]];
             if(iq.p2!=iqr.p2||iq.p4!=iqr.p4)
                printf("q2_list:momentum mode mismatch, %10d vs. %10d\n",
                     iq.iserial,iqr.iserial);
             int q_count=0;
             offset.push_back(list_new.size());
             for(int j=last_off;j<=i;j++)
             {
                if(list[j]->src.sub(list[j]->sink).iserial==iqr.iserial)
                {
                    list_new.push_back(new q2_mode(*list[j]));
                    q_count++;
                }
                delete list[j];
             }
             Npmode.push_back(q_count);
         }
         if(iseed==0)
         {
             Npmode.push_back(size_t);
             offset.push_back(last_off);
         }
         if(iseed<0)
         {
             int q_count=0;
             offset.push_back(list_new.size());
             for(int j=last_off;j<=i;j++)
             {
                 int is=momX.ind_2_mode(list[j]->sink);
                 if(list[j]->sink.iserial==momX.mode[is][sink_offset[is]].iserial)
                 {
                      list_new.push_back(new q2_mode(*list[j]));
                      q_count++;
                 }
                 delete list[j];
             }
             Npmode.push_back(q_count);
         }
         last_off=i+1;
      }
      Nq2=offset.size();
      if(iseed!=0)
      {
          list.resize(0);
          for(unsigned int i=0;i<list_new.size();i++)
              list.push_back(list_new[i]);
          list_new.resize(0);
      }
      size=list.size();
  }
  
  q2_list(std::vector<mom_mode> &src,std::vector<mom_mode> &sink,dimensions &desc,int cut_off=100,int iseed=0):desc(desc)
  { set(src,sink,cut_off,iseed); }
  q2_list(dimensions &desc):desc(desc)
  {
    mom_mode a(505050,1.0,desc);
    std::vector<mom_mode> list,list0;
    list0.push_back(a);list.push_back(a);
    set(list0,list);
  }
  
  void print(bool all=true){
  if(all==true)
   for(int i=0;i<Nq2;i++)
   for(int j=0;j<Npmode[i];j++)
   {
      printf("%4d%6d%6d:",i,j,offset[i]+j);
      list[offset[i]+j]->print(true);
   }
  else
  {
     for(int i=0;i<Nq2;i++)
     {
        printf("%4d%6d:",i,Npmode[i]);
        list[offset[i]]->print(true);
     }
  }  
  }

  ~q2_list()
  {
      for(unsigned int i=0;i<list.size();i++) delete list[i];
  }

  int find_q2(int ind_src,int ind_sink,double mass=1.0){
    mom_mode a(ind_src,mass,desc),b(ind_sink,mass,desc);
    q2_mode q2_tmp(a,b);
    return find_q2(q2_tmp.q2);  
  }
  int find_q(int ind_src,int ind_sink,double mass=1.0){
     mom_mode a(ind_src,mass,desc),b(ind_sink,mass,desc);
     q2_mode q2_tmp(a,b);
     int iq2=find_q2(q2_tmp.q2);
     for(int i=0;i<Npmode[iq2];i++)
       if(list[offset[iq2]+i]->src.iserial!=ind_src) continue;
       else
          if(list[offset[iq2]+i]->sink.iserial==ind_sink)return offset[iq2]+i;
      return -1;
  }
  
  int find_q2(double q2)
  {  return find_q2(q2,0,Nq2-1); }
  double p2i(int i)
  {  return list[offset[i]]->p2i(); }
  double p2f(int i)
  {  return list[offset[i]]->p2f(); }
  double q2(int i)
  {  return list[offset[i]]->q2; }
  
private:
  int find_q2(double q2,int ist,int ied){
  if(ied-ist<=1)
  {
    if(fabs(q2-list[offset[ist]]->q2)<1e-5) return ist;
    if(fabs(q2-list[offset[ied]]->q2)<1e-5) return ied;
    return -1;
  }
  double q2_tmp=list[offset[(ist+ied)/2]]->q2;
  if(fabs(q2-q2_tmp)<1e-5) return (ist+ied)/2;
  if(q2>q2_tmp) return find_q2(q2,(ist+ied)/2+1,ied);
  else return find_q2(q2,ist,(ist+ied)/2-1);  
  }
  
  static int comparer(q2_mode &a, q2_mode &b)
  {
       if(a.q2>b.q2) return 1;
       if(a.q2==b.q2)return 0;
       return -1;
  }
};

struct p2_list: q2_list
{
     p2_list(int max_mom,dimensions &desc,int iseed=0):q2_list(desc)
     {
         mom_set momX;momX.set();
         std::vector<mom_mode> list,list0;
         mom_mode Mtmp(0,0,0,1.0,desc);
         list0.push_back(Mtmp);
         for(int i=0;i<max_mom;i++)
         if(momX.mode[i].size()>0)
         {
            int imin=0,imax=momX.mode[i].size()-1;
            for(int j=imin;j<=imax;j++)
            {
               mom_mode mom_tmp(momX.mode[i][j],1.0,desc);
               list.push_back(mom_tmp);
            } 
         }
         set(list0,list,100,iseed);
     }
     // note that the random picked list just include one of the momentum for kinds of moemntum with the same p^2.
     
     int find_p2(momentum ind_sink)
     {  return find_q2(505050, ind_sink,1.0); }
     double p2(int i){ return p2f(i); }
};

}

#endif
