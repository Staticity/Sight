#include <utils/string_helpers.hpp>

#include <regex>

#include <codecvt>
#include <Windows.h>

namespace sight
{
    std::string string_cast(const std::wstring& wstr)
    {
        if (wstr.empty())
        {
            return {};
        }

        const int size_needed = WideCharToMultiByte(CP_UTF8, 0, &wstr[0], (int)wstr.size(), NULL, 0, NULL, NULL);
        std::string strTo(size_needed, 0);
        WideCharToMultiByte(CP_UTF8, 0, &wstr[0], (int)wstr.size(), &strTo[0], size_needed, NULL, NULL);
        return strTo;
    }

    std::wstring string_cast(const std::string& str)
    {
        if (str.empty())
        {
            return {};
        }
        
        const int size_needed = MultiByteToWideChar(CP_UTF8, 0, &str[0], (int)str.size(), NULL, 0);
        std::wstring wstrTo(size_needed, 0);
        MultiByteToWideChar(CP_UTF8, 0, &str[0], (int)str.size(), &wstrTo[0], size_needed);
        return wstrTo;
    }

    std::string duplicate_backslashes(const std::string& str)
    {
        return std::regex_replace(str, std::regex(R"(\\)"), R"(\\)");
    }
}
