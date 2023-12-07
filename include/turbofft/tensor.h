namespace turbofft{
namespace fft{
    template<
    typename DataType,
    size_t... Dimensions
    > class Tensor{
        private:
            DataType* data;
            size_t NumDimensions = sizeof...(Dimensions);
            size_t dimensions[sizeof...(Dimensions)]={Dimensions...};
        public:
            Tensor(){}

            size_t CalculateTotalElements(){
                return (Dimensions * ... * 1);
            }

            size_t CalculateTotalSize(){
                return (Dimensions * ... * 1) * sizeof(DataType);
            }
    };




}
}