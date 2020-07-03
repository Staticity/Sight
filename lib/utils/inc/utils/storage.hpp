#pragma once

namespace sight
{

    constexpr int STATIC_STORAGE = 0;
    constexpr int REF_STORAGE = 1;

    template <typename S, int N, int Options>
    struct Storage;
    
    template <typename S, int N>
    struct Storage<S, N, STATIC_STORAGE>
    {
        Storage()
        {}

        Storage(S* p)
        {
            std::copy(p, p + N, begin());
        }

        Storage(std::initializer_list<S> l)
        {
            std::copy(l.begin(), l.end(), begin());
        }

        template <int OtherOptions>
        Storage(const Storage<S, N, OtherOptions>& s)
        {
            std::copy(s.begin(), s.end(), begin());
        }
        
        S* begin() { return &v[0]; }
        const S* begin() const { return &v[0]; }

        S* end() { return begin() + N; }
        const S* end() const { return begin() + N; }

        S& operator[](int i) { return v[i]; }
        const S& operator[](int i) const { return v[i]; }

        S v[N];
    };
    
    template <typename S, int N>
    struct Storage<S, N, REF_STORAGE>
    {
        // Cannot implemented
        Storage();
        Storage(std::initializer_list<S> l);
        Storage(S* p)
        {
            v = p;
        }

        template <int OtherOptions>
        Storage(const Storage<S, N, OtherOptions>& s)
        {
            v = s.begin();
        }

        S* begin();
        S* end();
        S& operator[](int i);

        const S* begin() const { return v; }
        const S* end() const { return begin() + N; }
        const S& operator[](int i) const { return v[i]; }

        const S* v;
    };


} // namespace sight