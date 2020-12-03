#pragma once

#include <map>
#include <set>
#include <vector>

namespace sight
{
    class UnionFindDense
    {
    public:
        UnionFindDense(int capacity);

        void Find(int x);
        void Join(int a, int b);

        std::set<int> GroupIds() const;
        std::map<int, std::vector<int>> Groups() const;

    private:

        std::vector<int> m_rank;
        mutable std::vector<int> m_parent;
    };
}
