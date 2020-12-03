#pragma once

#include <string>

namespace sight
{
    std::wstring string_cast(const std::string& utf8Str);
    std::string string_cast(const std::wstring& utf16Str);
    std::string duplicate_backslashes(const std::string& str);
}
