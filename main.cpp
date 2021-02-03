#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include <memory>

using std::vector;
using std::cerr;
using std::ostream;
using Number = float;

namespace {
struct Array;
}


template <typename T>
static void destroyObject(T &object)
{
  object.~T();
}


namespace {
struct Value {
  enum class Type {
    number, character, array_ptr
  };

  Type type;

  union {
    Number number;
    char character;
    std::unique_ptr<Array> array_ptr;
  };

  Value()
  : type(Type::array_ptr), array_ptr(nullptr)
  {
  }

  Value(Number arg)
  : type(Type::number), number(arg)
  {
  }

  Value(char arg)
  : type(Type::character), character(arg)
  {
  }

  Value(const Value &arg)
  : type(arg.type)
  {
    createFrom(arg);
  }

  bool isNumber() const { return type == Type::number; }
  bool isCharacter() const { return type == Type::character; }

  Number asNumber() const
  {
    assert(type == Type::number);
    return number;
  }

  friend bool operator==(const Value &a, const Value &b)
  {
    if (a.type != b.type) {
      assert(false);
    }

    switch (a.type) {
      case Type::number:
        return a.number == b.number;
      case Type::character:
        assert(false);
      case Type::array_ptr:
        assert(false);
    }
    assert(false);
  }

  void destroy()
  {
    switch (type) {
      case Type::number:
        destroyObject(number);
        break;
      case Type::character:
        destroyObject(character);
        break;
      case Type::array_ptr:
        destroyObject(array_ptr);
        break;
    }
  }

  void createFrom(const Value &arg)
  {
    assert(type == arg.type);

    switch (type) {
      case Type::number:
        new (&number) auto(arg.number);
        break;
      case Type::character:
        new (&character) auto(arg.character);
        break;
      case Type::array_ptr:
        assert(false);
    }
  }

  Value& operator=(const Value &arg)
  {
    destroy();

    type = arg.type;

    createFrom(arg);

    return *this;
  }

  ~Value() { destroy(); }
};
}


namespace {
struct Array {
  vector<int> shape;
  vector<Value> values;

  friend bool operator==(const Array &a, const Array &b)
  {
    if (a.shape != b.shape) {
      if (a.values.size() == 1 && b.values.size() == 1) {
        return true;
      }
      assert(false);
    }

    if (a.values != b.values) {
      return false;
    }

    return true;
  }
};
}


static Array makeScalarArray(Value value)
{
  Array result;
  result.values = {value};
  return result;
}


static Array calculate(Number arg)
{
  return makeScalarArray(arg);
}


static Array calculate(const char *arg)
{
  Array result;
  int n = strlen(arg);
  result.shape = vector<int>{ n };
  result.values.resize(n);

  for (int i=0; i!=n; ++i) {
    result.values[i] = arg[i];
  }

  assert(result.values[0].isCharacter());
  return result;
}


static ostream& operator<<(ostream& , const Value &)
{
  assert(false);
}


template <typename T>
static ostream& operator<<(ostream& s, const vector<T> &v)
{
  s << "[";

  for (const auto &x : v) {
    s << " " << x;
  }

  s << " ]";
  return s;
}


static Array join(Array arg1, Array arg2)
{
  if (arg1.shape.size() == 0 && arg2.shape.size() == 0) {
    Array a;
    a.shape.resize(1);
    a.shape[0] = 2;
    a.values.resize(2);
    a.values[0] = arg1.values[0];
    a.values[1] = arg2.values[0];
    return a;
  }
  else if (arg1.shape.size() == 0 && arg2.shape.size() == 1) {
    Array a;
    a.shape.resize(1);
    a.shape[0] = arg2.shape[0] + 1;
    a.values.resize(a.shape[0]);
    a.values[0] = arg1.values[0];
    std::copy(arg2.values.begin(), arg2.values.end(), a.values.begin() + 1);
    return a;
  }
  else {
    cerr << "arg1.shape: " << arg1.shape << "\n";
    cerr << "arg2.shape: " << arg2.shape << "\n";
    assert(false);
  }
}


namespace {
struct Shape { };
struct Equal { };
struct Plus { };
struct Iota { };
struct Reduce { };
}


namespace {
template <typename T>
struct Primitive {
  T arg;
};
}


namespace {
template <typename T>
struct BoundRight {
  Array right;
};
}


static Array evaluate(BoundRight<Shape> arg)
{
  Array array = arg.right;
  Array result;
  result.shape = { int(array.shape.size()) };

  for (int x : array.shape) {
    result.values.push_back(Number(x));
  }

  return result;
}


static Array evaluate(BoundRight<Iota> arg)
{
  if (arg.right.shape.size() == 0) {
    int n = arg.right.values[0].asNumber();
    assert(n >= 0);
    Array result;
    result.shape = {n};

    for (int i=1; i<=n; ++i) {
      result.values.push_back(Number(i));
    }

    return result;
  }

  Array result;
  //result.shape = { int(
  assert(false);
}


template <typename T>
static BoundRight<T> join(Primitive<T>, Array array)
{
  return { array };
}


namespace {
template <typename T>
struct BoundBoth {
  Array left;
  Array right;
};
}


template <typename T>
static BoundBoth<T> join(Array left, BoundRight<T> equal_array)
{
  return { left, equal_array.right };
}


template <typename T>
static BoundBoth<T> join(Array arg1, BoundBoth<T> args)
{
  return { join(arg1, args.left), args.right };
}


template <typename T>
static Primitive<T> calculate(Primitive<T> (*)())
{
  return {};
}


template <typename T>
static Primitive<T> calculate(Primitive<T> arg)
{
  return arg;
}


static BoundRight<Equal> join(Primitive<Equal> left, BoundRight<Shape> right)
{
  return join(left, evaluate(right));
}


namespace {
template <typename Function, typename Operator>
struct BoundOperator {
};
}


static BoundRight<Reduce>
join(Primitive<Reduce> /*left*/, BoundRight<Iota> right)
{
  return {evaluate(right)};
}


static BoundRight<BoundOperator<Plus,Reduce>>
join(Primitive<Plus> /*left*/, BoundRight<Reduce> right)
{
  return {right.right};
}


template <typename Arg1, typename Arg2, typename ...Args>
static auto calculate(Arg1 arg1, Arg2 arg2, Args ...args)
{
  return join(calculate(arg1), calculate(arg2, args...));
}


template <typename Function>
static Array evaluateBinary(Array left, Array right, Function f)
{
  if (left.values.size() == 1) {
    if (right.values.size() == 1) {
      Array result;
      result.shape = vector<int>{3};

      result.values =
        vector<Value>{ Value(f(left.values[0], right.values[0])) };

      return result;
    }
    assert(false);
  }

  if (left.shape == right.shape) {
    Array result;
    result.shape = left.shape;
    result.values.resize(left.values.size());
    auto in1 = left.values.begin();
    auto in2 = right.values.begin();
    auto out = result.values.begin();

    while (out != result.values.end()) {
      *out++ = f(*in1++, *in2++);
    }

    return result;
  }

  cerr << "left.values: " << left.values << "\n";
  assert(false);
}


static Array evaluate(BoundBoth<Equal> arg)
{
  auto f = [](auto a, auto b){ return Number(a == b); };
  return evaluateBinary(arg.left, arg.right, f);
}


static Array evaluate(BoundBoth<Plus> arg)
{
  auto f = [](const Value& a, const Value& b)
  {
    if (a.isNumber() && b.isNumber()) {
      return Value(a.asNumber() +b.asNumber());
    }
    else {
      assert(false);
    }

    return Value();
  };

  return evaluateBinary(arg.left, arg.right, f);
}


static Array evaluate(Array array)
{
  return array;
}


static Array evaluate(BoundRight<BoundOperator<Plus,Reduce>> arg)
{
  const Array &right = arg.right;

  if (right.shape.size() == 1) {
    Number total = 0;

    for (auto &x : right.values) {
      if (!x.isNumber()) {
        assert(false);
      }

      total += x.asNumber();
    }

    return makeScalarArray(total);
  }

  assert(false);
}


namespace {
struct Placeholder {
  template <typename ...Args>
  auto operator()(Args ...args)
  {
    return evaluate(calculate(args...));
  }

  static Primitive<Shape>  shape()  { return {}; }
  static Primitive<Equal>  equal()  { return {}; }
  static Primitive<Plus>   plus()   { return {}; }
  static Primitive<Iota>   iota()   { return {}; }
  static Primitive<Reduce> reduce() { return {}; }
};
}


int main()
{
  Placeholder _;
  assert(_(1) == _(1));
  assert(!(_(2) == _(3)));
  assert(_(1,2) == _(1,2));
  assert(_(1,2,3) == _(1,2,3));
  assert(_(1,2,3).shape == vector<int>{3});
  assert(_(_.shape, 1,2,3) == _(3));
  assert(_(3, _.equal, _.shape, 1,2,3) == _(1));
  assert(_(_.shape,1,2,3).shape == vector<int>{1});
  assert(_(1,2,3, _.equal, 1,2,3) == _(1,1,1));
  assert(_(2, _.plus, 2) == _(4));
  assert(_(_.shape, "hello") == _(5));
  assert(_(_.iota, 3) == _(1,2,3));
  assert(_(_.plus, _.reduce, _.iota, 3) == _(6));
}
