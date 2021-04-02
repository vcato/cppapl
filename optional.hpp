#include <utility>
#include <cassert>
#include <new>
#include "lifetime.hpp"


template <typename T>
class Optional {
  public:
    Optional(T arg)
    : has_value(true),
      value(std::move(arg))
    {
    }

    Optional(Optional &&arg)
    : has_value(arg.has_value)
    {
      if (has_value) {
        createObject(value, std::move(arg.value));
      }
    }

    Optional()
    : has_value(false)
    {
    }

    ~Optional()
    {
      if (has_value) {
        value.~T();
      }
    }

    bool operator!() const
    {
      return !has_value;
    }

    explicit operator bool() const
    {
      return has_value;
    }

    Optional& operator=(const T &arg)
    {
      if (has_value) {
        value = arg;
      }
      else {
        createObject(value, arg);
        has_value = true;
      }

      return *this;
    }

    T &operator*()
    {
      assert(has_value);
      return value;
    }

    T *operator->()
    {
      assert(has_value);
      return &value;
    }

    friend bool operator==(const Optional &a, const Optional &b)
    {
      if (a.has_value != b.has_value) {
        return false;
      }

      if (!a.has_value) {
        return true;
      }

      assert(b.has_value);
      return a.value == b.value;
    }

    friend bool operator!=(const Optional &a, const Optional &b)
    {
      return !operator==(a,b);
    }

  private:
    bool has_value;

    union {
      T value;
    };
};
