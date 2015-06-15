#include <iterator> /* std::iterator_traits */


namespace utils {


    template <class InputIterator, class OutputIterator>
    void exclusivePrefixSum(InputIterator first, InputIterator last, OutputIterator result){

	typedef typename std::iterator_traits<InputIterator>::value_type IterType;

	IterType value = 0;

	while(first != last){
	    IterType prevValue = value;
	    value   = value + *first;
	    *result = prevValue;
	    result++; first++;
	}
  
    }

} /* utils */
