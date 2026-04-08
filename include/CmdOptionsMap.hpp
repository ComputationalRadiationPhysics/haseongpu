/**
 * Copyright 2026 Tim Hanel
 */
#pragma once

#include <logging.hpp>
#include <types.hpp>

#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
namespace fs = std::filesystem;

std::string_view requireValue(int& i, int argc, char** argv, std::string_view opt);

template<typename T>
T sToX(std::string const& s);

template<>
inline int sToX<int>(std::string const& s)
{
    return std::stoi(s);
}

template<>
inline double sToX<double>(std::string const& s)
{
    return std::stod(s);
}

template<>
inline bool sToX<bool>(std::string const& s)
{
    if(s == "true" || s == "1")
        return true;
    if(s == "false" || s == "0")
        return false;
    throw std::runtime_error("invalid bool");
}

template<>
inline std::string sToX<std::string>(std::string const& s)
{
    return s;
}

template<>
inline fs::path sToX<fs::path>(std::string const& s)
{
    return fs::path{s};
}

template<>
inline unsigned sToX<unsigned>(std::string const& s)
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
            if(hasValue)
            {
                return sToX<T>(value);
            }
            return sToX<T>(defaultValue);
        }

        template<typename T>
        T as()
        {
            if(hasValue)
            {
                return sToX<T>(value);
            }
            return sToX<T>(defaultValue);
        }
    };

    std::unordered_map<std::string, Option> options;

    auto& operator[](std::string name)
    {
        return options[name];
    }

    CmdOptionsMap();

    void add(std::string name, std::string description, std::string defaultValue);

    bool isExplicitlySet(std::string const& name) const;
    bool contains(std::string const& name) const;

    template<typename T>
    T as(std::string const& name) const
    {
        auto it = options.find(name);
        if(it == options.end())
        {
            throw std::runtime_error("Unknown option: " + name);
        }
        return it->second.as<T>();
    }

    void set(std::string const& name, std::string value);

    void printDescription() const;
    void checkRequired();

private:
    void construct();
};
