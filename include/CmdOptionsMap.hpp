/**
* Copyright 2026 Tim Hanel
*/
#pragma once

#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <logging.hpp>
#include <types.hpp>
namespace fs = std::filesystem;

std::string_view requireValue(int& i, int argc, char** argv, std::string_view opt);

template<typename T>
T sToX(const std::string& s);

template<>
inline int sToX<int>(const std::string& s)
{
    return std::stoi(s);
}

template<>
inline double sToX<double>(const std::string& s)
{
    return std::stod(s);
}

template<>
inline bool sToX<bool>(const std::string& s)
{
    if (s == "true" || s == "1") return true;
    if (s == "false" || s == "0") return false;
    throw std::runtime_error("invalid bool");
}

template<>
inline std::string sToX<std::string>(const std::string& s)
{
    return s;
}

template<>
inline fs::path sToX<fs::path>(const std::string& s)
{
    return fs::path{s};
}

template<>
inline unsigned sToX<unsigned>(const std::string& s)
{
    return static_cast<unsigned>(std::stoul(s));
}

struct CmdOptionsMap
{
    struct Option
    {
        std::string name;
        std::string description;
        std::string value;
        std::string defaultValue;
        bool hasValue = false;

        void set(std::string newValue);

        template<typename T>
        T as() const
        {
            if (hasValue) {
                return sToX<T>(value);
            }
            return sToX<T>(defaultValue);
        }

        template<typename T>
        T as()
        {
            if (hasValue) {
                return sToX<T>(value);
            }
            return sToX<T>(defaultValue);
        }
    };

    std::unordered_map<std::string, Option> options;
    auto& operator[](std::string name){
        return options[name];
    }
    CmdOptionsMap();

    void add(
        std::string name,
        std::string description,
        std::string defaultValue);

    bool isExplicitlySet(const std::string& name) const;
    bool contains(const std::string& name) const;

    template<typename T>
    T as(const std::string& name) const
    {
        auto it = options.find(name);
        if (it == options.end()) {
            throw std::runtime_error("Unknown option: " + name);
        }
        return it->second.as<T>();
    }

    void set(const std::string& name, std::string value);

    void printDescription() const;
    void checkRequired();

private:
    void construct();
};