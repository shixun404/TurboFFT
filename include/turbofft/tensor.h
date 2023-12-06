namespace turbofft{
namespace fft{
    template<
    typename DataType,
    int Dim,
    int *Dims
    > class Tensor{
        private:
            DataType* data;
            int dim = Dim;
            int dims[Dim]; //initiliaze dims with dim
        public:
            Tensor(int dims*){
                #pragma unroll
                for(int i = 0; i < Dim; ++i){
                    dims[i] = Dims[i];
                }
            }
    }


}
}