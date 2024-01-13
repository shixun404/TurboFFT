namespace turbofft{
    template<
    typename DataType,
    size_t... Dimensions
    > class Tensor{
        private:
            DataType* data;
            // size_t NumDimensions = sizeof...(Dimensions);
            size_t NumDimensions = 1;
            size_t dimensions[1]={1};
            // size_t stride[sizeof...(Dimensions)]={Dimensions...};
        public:
            Tensor(DataType* data){
                this->data = data;
            }

            size_t CalculateTotalElements(){
                return 1;
                // return (Dimensions * ... * 1);
            }

            size_t CalculateTotalSize(){
                return 1;
                // return (Dimensions * ... * 1) * sizeof(DataType);
            }
    };
    // template<typename DataType>
    // class Tensor<DataType, 1, 2>{
    //     private:
    //         DataType* data;
    //         size_t NumDimensions = 2;
    //         size_t dimensions[2]={1, 2};
    //         // size_t stride[sizeof...(Dimensions)]={Dimensions...};
    //     public:
    //         Tensor(DataType* data){
    //             this->data = data;
    //         }

    //         size_t CalculateTotalElements(){
    //             return (1 * 2);
    //         }

    //         size_t CalculateTotalSize(){
    //             return (1 * 2) * sizeof(DataType);
    //         }
    // };
}