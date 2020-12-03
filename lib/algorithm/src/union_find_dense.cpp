#include <algorithm/union_find.hpp>

#include <stack>
#include <numeric>

namespace sight
{
    UnionFindDense::UnionFindDense(int capacity)
        : m_parent(capacity)
        , m_rank(capacity, 0)
    {
        std::iota(m_parent.begin(), m_parent.end(), 0);
    }

    int UnionFindDense::Find(int x)
    {
        // Traverse up the parent links to
        // find the representative.
        //
        // Store a history of the child nodes visited.
        stack<int> s;
        while (m_parent[x] != x)
        {
            s.push(x);
            x = m_parent[x];
        }

        // Store the max depth as the rank of this node.
        m_rank[x] = std::max(m_rank[x], s.size());

        // Label all of the child nodes to the discovered parent.
        while (!s.empty())
        {
            m_parent[s.top()] = x;
            s.pop();
        }

        return x;
    }

    void UnionFindDense::Join(int a, int b)
    {
        const int pa = Find(a);
        const int pb = Find(b);

        // We attach to the node which has the
        // larger depth. This keeps the depth
        // at a small size.
        if (m_rank[pa] < m_rank[pb])
        {
            m_parent[pa] = pb;
        }
        else
        {
            m_parent[pb] = pa;
        }
    }

    std::set<int> UnionFindDense::GroupIds() const
    {
        std::set<int> ids;
        for (int x : m_parent)
        {
            ids.insert(x):
        }
        return ids;
    }

    std::map<int, std::vector<int>> Groups() const
    {
        std::map<int, std::vector<int>> groups;
        for (int i = 0; i < m_parent.size(); ++i)
        {
            groups[i].push_back(Find(i));
        }
        return groups;
    }

}
