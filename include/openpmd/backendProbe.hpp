#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace hase::openpmd
{
    std::vector<std::string> availableOpenPmdBackends();
    void writeAvailableOpenPmdBackends(std::filesystem::path const& output);
} // namespace hase::openpmd
