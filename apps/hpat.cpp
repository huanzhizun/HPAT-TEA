/*************************************************************************
    > File Name: random.cpp
    > Author: huanzhizun
    > Mail: huznzhiun888@126.com 
    > Created Time: Fri 09 Mar 2018 03:32:33 PM CST
 ************************************************************************/
#include "src/core.hpp"
#include<queue>
using namespace std;
const double eps=1e-9;
const int logd=15;
const int logdd=20;
//typedef pair< pair<int,char>,long long> ran_type;
typedef pair<char ,int> bin_type;
typedef pair<pair<int,long long>,int> task_type;
bin_type* per_bin;
int* bin_begin;
int concurrency=16;
task_type *task;
int *be,*en,*en_t;
char *bin_num;
long long *bin;
vt *weight;
int *alias_i;
double *alias_d;
long long *ver_bin;
StdRandNumGenerator* rand_gen;
//vector<int> alias_i;
//vector<double> alias_d;
//vector<ran_type>ver_bin,ver_tmp;
//vector<ran_type>ver_tmp;
int vertices,degree;
long long edges;
long long vv[32];
int rap[32];
long long base_l=0,base_v=0;
/*struct TimeEdge {
	ut u, v;
	ut s;
};
TimeEdge *edge;
*/
ut* edge;
char *mmap_index;
char *mmap_edges;
long long mmap_edges_size;
void init_edge(string edge_file_name){
	edge=(ut *)malloc((long long)3*sizeof(ut)*(edges+10));
	double start_time = get_time();
	long bytes = file_size(edge_file_name);
	int fd = open(edge_file_name.c_str(), O_RDONLY);
	posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
	char * buffers = (char *)mmap(0, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	long read_bytes = 0;
	while(true) {
		long r_bytes = read(fd, buffers + read_bytes, IOBUFFER);
		if (r_bytes <= 0) break;
		read_bytes += r_bytes;
	}
	close(fd);
	assert(read_bytes == bytes);
    for(long long i=0;i<(long long)3*edges;i++){
        ut *p=(ut *)(buffers+i*sizeof(ut));
        edge[i]=*p;
    }
}
void init_vertice(string input_vertice){
	long bytes = file_size(input_vertice);
	int fd = open(input_vertice.c_str(), O_RDONLY);
	posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
	char * buffers = (char *)mmap(0, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	double start_time = get_time();
	long read_bytes = 0;
	while(true) {
		long r_bytes = read(fd, buffers + read_bytes, IOBUFFER);
		if (r_bytes <= 0) break;
		read_bytes += r_bytes;
	}
	close(fd);
	assert(read_bytes == bytes);
    weight=(vt *)malloc(sizeof(vt)*(edges+10));
    bin_num=(char *)malloc(sizeof(char)*(degree+10));
    bin_begin=(int *)malloc(sizeof(int)*(degree+10));
    per_bin=(bin_type*)malloc(sizeof(bin_type)*(degree*logd+10)); 
    task=(task_type*)malloc(sizeof(task_type)*(edges*3+10)); 
    be=(int *)malloc(sizeof(int)*(vertices+10));
    en=(int *)malloc(sizeof(int)*(vertices+10));
    bin=(long long *)malloc(sizeof(long long)*(edges+10));
    alias_d=(double *)malloc(sizeof(double)*((long long)edges*logdd+10));
    alias_i=(int *)malloc(sizeof(int)*((long long)edges*logdd+10));
    ver_bin=(long long *)malloc(sizeof(long long)*((long long)edges*logd+10));
    for(ut i=0;i<vertices;i++){
        ut *p=(ut *)(buffers+i*2*sizeof(int));
        be[i+1]=*p;
        p=(ut *)(buffers+(i*2+1)*sizeof(int));
        en[i+1]=*p;
    }
 /*	int index_fd = open(input_vertice.c_str(), O_RDONLY);
    mmap_index = (char *)mmap(0, sizeof(ut) * vertices, PROT_WRITE, MAP_PRIVATE, index_fd, 0);
    label = (ut *)mmap_index;
    close(index_fd);
    madvise(mmap_index, sizeof(ut) *vertices, MADV_WILLNEED)
*/
}
void init(){
	weight[0]=0;
    for(int i=0;i<edges;i++){
        weight[i+1]=edge[(long long)3*i+1];
        weight[i+1]+=weight[i];
    }
    int all=0;
    for(int i=1;i<=degree+2;i++){
        bin_begin[i]=all;
        int sum=0;
        for(int j=30;j>=0;j--){
            if((i>>j)&1){
                per_bin[all++]=make_pair(j,sum);
                sum+=(1<<j);
            }
        }
        bin_num[i]=all-bin_begin[i];
    }
    for(int i=0;i<32;i++){
        rap[i]=1<<i;
    }
}
void alias(long long base_l,int base_lt,int len){
    queue<int>smaller;
    queue<int>larger;
    for(int i=0;i<len;i++){
        if (alias_d[base_l+i] < 1.0) smaller.push(i);
        else larger.push(i);
    }
    while(!smaller.empty()&&!larger.empty()){
        int s_t=smaller.front();
        smaller.pop();
        int l_t=larger.front();
        larger.pop();   
        alias_i[base_l+s_t]=base_lt+l_t;
        alias_d[base_l+l_t]-=1-alias_d[base_l+s_t];
        if(alias_d[l_t]<1){
            smaller.push(l_t);
        }
        else{
            larger.push(l_t);
        }
    }
    while(!smaller.empty()){
        int w=smaller.front();
        smaller.pop();
        alias_d[base_l+w]=1;
    }
    while(!larger.empty()){
        int w=larger.front();
        larger.pop();
        alias_d[base_l+w]=1;
    }
}
void update_task(task_type num_tmp){
    int w=num_tmp.first.first;
    long long alias_base=num_tmp.first.second;
    int be_base=num_tmp.second;
    double tmp=0;
    for(int u=0;u<w;u++){
        tmp+=edge[(long long)3*(be_base+u)+1];
    }
    for(int u=0;u<w;u++){
        alias_d[alias_base+u]=edge[(long long)3*(be_base+u)+1]/tmp*w;
        alias_i[alias_base+u]=be_base+u;
    }
    alias(alias_base,be_base,w);
}
/*
void update_alias(int q,int len,int ver){
    Queue<task_type> tasks(concurrency);
    std::vector<std::thread> threads;
//    #pragma omp parallel for num_threads(concurrency)
    for (int i = 0; i < concurrency; i++){
        threads.emplace_back([&](){
        while (true) {
            task_type num_tmp = tasks.pop();
            if(num_tmp.second==-1) break;
            update_task(num_tmp,ver);
            }
        });
    }
    for(int j=0;j<q;j++){
        int w=1<<j;
        for(int k=0;k<=len-w;k+=w){
            tasks.push(make_pair(j,k));
                    //alias(base_l,be[i]+k,w);
                    //base_l+=w;
        }
    }
    for (int i = 0; i < concurrency; i++){
       tasks.push(make_pair(-1,-1));
    }
    for (int i = 0; i < concurrency; i++) threads[i].join();
}
void update_alias_tmp(int q,int len,int ver){
    int task_num=0;
    for(int j=q-1;j>=0;j--){
        int w=rap[j];
        for(int k=0;k<=len-w;k+=w){
            task[task_num++]=make_pair(j,k);
        }
    }
    #pragma omp parallel for schedule(guided,concurrency)
    for(int i=0;i<task_num;i++){
        update_task(task[i],ver);
    }
}*/
void update_bin(int ver){
    for(int j=be[ver];j<=en[ver];j++){
        int v=j-be[ver]+1;
        bin[j]=base_v;
        int w=bin_num[v];
        //#pragma omp parallel for schedule(dynamic,concurrency)
//        #pragma omp parallel for num_threads(concurrency)
        for(int k=0;k<w;k++){
            ver_bin[base_v+k]=vv[per_bin[bin_begin[v]+k].first];
        }
        base_v+=w;
    }
}
void init_alias(){
    double no_time=0;
    int task_num=0;
    for(int i=1;i<=vertices;i++){
//        if(i>=vertices_o-10)    
        if(be[i]!=-1){
            int len=en[i]-be[i]+1;
            long long w=1;
            int q=0;
            for(q=0;w<=len;q++){
                vv[q]=base_l;
                base_l+=len-(len%w);
                w=w<<1;
            }
            for(int j=0;j<q;j++){
                int w=rap[j];
                for(int k=0;k<=len-w;k+=w){
                    task[task_num++]=make_pair(make_pair(w,vv[j]+k),be[i]+k);
                }
            }
  //          update_alias_tmp(q,len,i);
            update_bin(i);
//            #pragma omp parallel for schedule(dynamic,concurrency)
            /*
            for(int j=be[i];j<=en[i];j++){
                int v=j-be[i]+1;
                int w=0;
                bin[j]=base_v;
                bin_num[j]=0;
                for(int k=31;k>=0;k--){
                    if((v>>k)&1){
                        bin_num[j]++;
                        //ver_bin.push_back(make_pair(make_pair(w+be[i],k),vv[k]+w));
                        ver_bin[base_v]=make_pair(make_pair(w+be[i],k),vv[k]+w);
                        w+=1<<k;
                        base_v++;
                    }
                }
            }*/
        }
    }
    #pragma omp parallel for schedule(dynamic,concurrency)
    for(int i=0;i<task_num;i++){
        update_task(task[i]);
    }
}
void update(int u,int v,int len,int thread_num){
    if(u==-1||v==-1||len==0){
        return ;
    }
    if(be[u]==-1){
        return ;
    }
    vt tmp=weight[v+1];
    //vt sum=tmp;
    vt sum=weight[be[u]];
    long long ran=rand_gen[thread_num].gen(tmp-sum)+1+sum;
    int le=0, ri=bin_num[v-be[u]+1]-1;
    int mid;
    while(le<=ri){
        mid=(le+ri)/2;
        int w=be[u]+per_bin[mid+bin_begin[v-be[u]+1]].second;
        //long long ran_t=weight[w.first.first];
        if(weight[w]<ran) le=mid+1;
        else ri=mid-1;
    }
    bin_type now_bin=per_bin[ri+bin_begin[v-be[u]+1]];
    long long alias_pos=ver_bin[bin[v]+ri]+now_bin.second;
    int begin_pos=be[u]+now_bin.second;
    double now=(double)(ran-weight[begin_pos])/(weight[begin_pos+rap[now_bin.first]]-weight[begin_pos])*rap[now_bin.first];
    int x=(int)(now-eps);
    double y=now-x;
    long long trans=begin_pos+x;
    if(y>alias_d[alias_pos+x]+eps) trans=alias_i[alias_pos+x];
    update(edge[trans*3],edge[3*trans+2],len-1,thread_num);
}
int main(int argc,char **argv)
{
	cmdline::parser cmd;
	cmd.add<std::string>("input", 'i', "input file", true, "");
	cmd.add<int>("number", 'n', "number of walkers", true, 0);
	cmd.add<int>("length", 'l', "length of walkers", true, 0);

	cmd.parse_check(argc, argv);
	string input = cmd.get<string>("input");
    ut num = cmd.get<int>("number");
    ut len = cmd.get<int>("length");

	printf("Begin \n");
	string input_edge=input+"-edge.info";
	string input_vertice=input+"-vertice.info";
	string input_info=input+".info";
	int fin = open(input_info.c_str(), O_RDONLY);
    int tmp_edges;
	if(read(fin, (char*)&vertices, sizeof(ut)) != sizeof(ut)) {
		printf("Format error\n");
		exit(-1);
	}
	if(read(fin, (char*)&tmp_edges, sizeof(ut)) != sizeof(ut)) {
		printf("Format error\n");
		exit(-1);
	}
    edges=tmp_edges;
	if(read(fin, (char*)&degree, sizeof(ut)) != sizeof(ut)) {
		printf("Format error\n");
		exit(-1);
	}
	printf("%d %d %d\n",vertices,edges,degree);
	init_vertice(input_vertice);
	init_edge(input_edge);
    init();
    srand((unsigned)time(NULL));
    init_alias();
    rand_gen=new StdRandNumGenerator[concurrency];
    {
        #pragma omp parallel for schedule(guided,concurrency)
        for(int i=1;i<=num;i++){
            int thread_num=omp_get_thread_num();
            long long w=rand_gen[thread_num].gen(edges);
             update(edge[w*3],edge[w*3+2],len,thread_num);
        }   
    }
    munmap(edge, mmap_edges_size);
	return 0;
}
