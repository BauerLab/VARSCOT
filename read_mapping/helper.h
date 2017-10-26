#include <type_traits>

using namespace std;
using namespace seqan;

template <typename T>
inline void printVector(vector<T> const & v, ostream & stream = cout)
{
    stream << "[";
    for (unsigned i = 0; i < v.size(); ++i)
    {
        stream << (is_same<T, uint8_t>::value ? (unsigned) v[i] : v[i])
               << ((i < v.size() - 1) ? ", " : "");
    }
    stream << "]" << endl;
}

template <typename TText, typename TRand>
void generateText(TText & text, TRand & rng, unsigned length)
{
    typedef typename Value<TText>::Type TChar;

    int minChar = MinValue<TChar>::VALUE;
    unsigned alphabetSize = ValueSize<TChar>::VALUE;

    resize(text, length);

    for (unsigned i = 0; i < length; ++i)
        text[i] = rng() % alphabetSize - minChar;
}

// TODO: the same as trivialSearch but only takes a unidirectional iter. both methods should be merged to one!
template <typename TDelegate, typename TText, typename TIndex,
        typename TIndexSpec, typename TText2, typename TNeedleIter>
inline void trivialUniSearch(TDelegate & delegate,
                          Iter<Index<TText, TIndex>, VSTree<TopDown<TIndexSpec> > > it,
                          TText2 const & needle,
                          TNeedleIter const needleIt,
                          uint8_t const errorsLeft,
                          bool const indels
#ifdef ENABLE_DEBUG_MACRO
        , std::string str = ""
#endif
)
{
    if (errorsLeft == 0)
    {
        if (goDown(it, suffix(needle, position(needleIt, needle))))
        {
            delegate(it
#ifdef ENABLE_DEBUG_MACRO
                    , str + toCString((CharString) suffix(needle, position(needleIt, needle)))
#endif
            );
        }
    }
    else
    {
        if (atEnd(needleIt, needle))
        {
            delegate(it
#ifdef ENABLE_DEBUG_MACRO
                    , str
#endif
            );
            if (!(indels && errorsLeft > 0))
                return;
        }
        if (indels && !atEnd(needleIt, needle))
        {
            // Insertion
            trivialUniSearch(delegate, it, needle, needleIt + 1, errorsLeft - 1, indels
#ifdef ENABLE_DEBUG_MACRO
                    , str + "+"
#endif
            );
        }
        if (goDown(it))
        {
            do
            {
                // Match / Mismatch
                if (!atEnd(needleIt, needle))
                {
                    uint8_t delta = !ordEqual(parentEdgeLabel(it), value(needleIt));
                    trivialUniSearch(delegate, it, needle, needleIt + 1, errorsLeft - delta, indels
#ifdef ENABLE_DEBUG_MACRO
                            , str + toCString((CharString) parentEdgeLabel(it))
#endif
                    );
                }

                if (indels)
                {
                    // Deletion
                    trivialUniSearch(delegate, it, needle, needleIt, errorsLeft - 1, indels
#ifdef ENABLE_DEBUG_MACRO
                            , str + "-"
#endif
                    );
                }
            } while (goRight(it));
        }
    }
}

// TODO: move to unidirectional indices for better performance comparison!
template <typename TDelegate, typename TText, typename TIndex,
          typename TIndexSpec, typename TText2, typename TNeedleIter>
inline void trivialSearch(TDelegate & delegate,
                          Iter<Index<TText, TIndex>, VSTree<TopDown<TIndexSpec> > > it,
                          TText2 const & needle,
                          TNeedleIter const needleIt,
                          uint8_t const errorsLeft,
                          bool const indels
                          #ifdef ENABLE_DEBUG_MACRO
                              , std::string str = ""
                          #endif
)
{
    if (errorsLeft == 0)
    {
        if (goDown(it, suffix(needle, position(needleIt, needle)), Rev()))
        {
            delegate(it
                #ifdef ENABLE_DEBUG_MACRO
                    , str + toCString((CharString) suffix(needle, position(needleIt, needle)))
                #endif
            );
        }
    }
    else
    {
        if (atEnd(needleIt, needle))
        {
            delegate(it
                #ifdef ENABLE_DEBUG_MACRO
                    , str
                #endif
            );
            if (!(indels && errorsLeft > 0))
                return;
        }
        if (indels && !atEnd(needleIt, needle))
        {
            // Insertion
            trivialSearch(delegate, it, needle, needleIt + 1, errorsLeft - 1, indels
                #ifdef ENABLE_DEBUG_MACRO
                      , str + "+"
                #endif
            );
        }
        if (goDown(it, Rev()))
        {
            do
            {
                // Match / Mismatch
                if (!atEnd(needleIt, needle))
                {
                    uint8_t delta = !ordEqual(parentEdgeLabel(it, Rev()), value(needleIt));
                    trivialSearch(delegate, it, needle, needleIt + 1, errorsLeft - delta, indels
                        #ifdef ENABLE_DEBUG_MACRO
                            , str + toCString((CharString) parentEdgeLabel(it))
                        #endif
                    );
                }

                if (indels)
                {
                    // Deletion
                    trivialSearch(delegate, it, needle, needleIt, errorsLeft - 1, indels
                        #ifdef ENABLE_DEBUG_MACRO
                            , str + "-"
                        #endif
                    );
                }
            } while (goRight(it, Rev()));
        }
    }
}

template <typename TChar>
inline void _countTrivialSearch(unsigned const needleLength,
                                unsigned const needlePos,
                                uint8_t const errorsLeft, bool const indels,
                                unsigned long long & counts)
{
    if (errorsLeft == 0 || needleLength == needlePos)
    {
        counts += needleLength - needlePos;
        return;
    }

    // Match / Mismatch
    unsigned long long countsMatch = 1, countsMismatchOrInsertion = 0;
    _countTrivialSearch<TChar>(needleLength, needlePos + 1, errorsLeft, indels, countsMatch);
    _countTrivialSearch<TChar>(needleLength, needlePos + 1, errorsLeft - 1, indels, countsMismatchOrInsertion);
    counts += countsMatch + (ValueSize<TChar>::VALUE - 1) * (countsMismatchOrInsertion + 1);

    if (indels)
    {
        // Deletion
        unsigned long long countsDeletion = 0; // TODO: do we count the deletion itself?
        _countTrivialSearch<TChar>(needleLength, needlePos, errorsLeft - 1, indels, countsDeletion);
        counts += ValueSize<TChar>::VALUE * countsDeletion;
        // Insertion
        counts += countsMismatchOrInsertion; // TODO: do we count the deletion itself?
    }
}

template <typename TChar>
inline unsigned long long countTrivialSearch(unsigned int needleLength, uint8_t const errors, bool const indels)
{
    unsigned long long counts = 0;
    _countTrivialSearch<TChar>(needleLength, 0, errors, indels, counts);
    return counts;
}
