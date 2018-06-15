#include <type_traits>

#include <seqan/index.h>

using namespace seqan;

namespace seqan {

template <typename TChar, typename TOwner>
struct SAValue<StringSet<String<TChar>, TOwner > >
{
    typedef Pair<uint16_t, uint32_t> Type;
};

template <typename TChar, typename TOwner>
struct SAValue<String<TChar, TOwner > >
{
    typedef uint32_t Type;
};

};

// Index type
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
typedef BidirectionalIndex<FMIndex<void, TMyFastConfig> > TIndexConfig;
