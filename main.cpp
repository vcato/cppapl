#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include <memory>

#define ADD_EACH 0
#define CHANGE_JOINING 0

#define CHANGE_BOUND_RIGHT 0
  // I think BoundRight can have a vector<Value> for the left

#define CHANGE_EVALUATE_NUMBER 0
  // Evaluating a value should yield a value.

using std::vector;
using std::cerr;
using std::ostream;
using Number = float;

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

  Value(const Array &arg)
  : type(Type::array_ptr), array_ptr(std::make_unique<Array>(arg))
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
        assert(arg.array_ptr);
        new (&array_ptr) auto(std::make_unique<Array>(*arg.array_ptr));
        break;
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


static ostream& operator<<(ostream& , const Value &)
{
  assert(false);
}


namespace {
struct Array {
  vector<int> shape;
  vector<Value> values;

  friend bool operator==(const Array &a, const Array &b)
  {
    if (a.shape != b.shape) {
      if (a.values.size() == 1 && b.values.size() == 1) {
        return a.values[0] == b.values[0];
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


#if !CHANGE_BOUND_RIGHT
namespace {
template <typename T>
struct BoundBoth {
  Array left;
  Array right;
};
}
#endif


namespace {
template <typename T>
struct Operator {
};
}


namespace {
template <typename Function, typename Operator>
struct BoundOperator {
};
}


static Array makeScalarArray(Value value)
{
  Array result;
  result.values = {value};
  return result;
}


static Array makeCharArray(const char *arg)
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


namespace {
template <typename T>
struct BoundRight {
  Array right;
#if CHANGE_BOUND_RIGHT
  vector<Value> left;
#endif
};
}


namespace {
template <typename T>
struct Primitive {
};
}


template <typename Left, typename Right>
struct Atop {
};


namespace {
struct Shape;
struct First;
#if ADD_EACH
struct Each;
#endif
struct Equal;
struct Plus;
struct Times;
struct Iota;
struct Reduce;
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


#if CHANGE_EVALUATE_NUMBER
static Value evaluate(Number arg)
{
  return arg;
}


static Value evaluate(Value arg)
{
  return arg;
}
#else
static Array evaluate(Number arg)
{
  return makeScalarArray(arg);
}
#endif


#if CHANGE_JOINING
static Array evaluate(vector<Value> arg)
{
  Array result;
  result.shape = { int(arg.size()) };
  result.values = arg;
  return result;
}
#endif


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
  assert(false);
}


static Array evaluate(BoundRight<First> arg)
{
  if (arg.right.shape.size() == 0) {
    assert(false);
  }

  if (arg.right.shape.size() == 1) {
    return makeScalarArray(arg.right.values[0]);
  }

  assert(false);
}


template <typename A, typename B>
static Atop<A,B> evaluate(Atop<A,B> arg)
{
  return arg;
}


template <typename T>
static Primitive<T> evaluate(Primitive<T> (*)())
{
  return {};
}


static Operator<Reduce> evaluate(Primitive<Reduce> (*)())
{
  return {};
}


#if ADD_EACH
static Operator<Each> evaluate(Primitive<Each> (*)())
{
  return {};
}
#endif


#if !CHANGE_BOUND_RIGHT
static Array evaluate(BoundBoth<Equal> arg)
{
  auto f = [](auto a, auto b){ return Number(a == b); };
  return evaluateBinary(arg.left, arg.right, f);
}
#endif


#if !CHANGE_BOUND_RIGHT
template <typename T, typename G>
static Array evaluateNumberBinary(BoundBoth<T> arg, const G &g)
{
  auto f = [&](const Value& a, const Value& b)
  {
    if (a.isNumber() && b.isNumber()) {
      return Value(g(a.asNumber(), b.asNumber()));
    }
    else {
      assert(false);
    }

    return Value();
  };

  return evaluateBinary(arg.left, arg.right, f);
}
#endif


#if !CHANGE_BOUND_RIGHT
static Array evaluate(BoundBoth<Plus> arg)
{
  return evaluateNumberBinary(arg, [](Number a, Number b) { return a + b; });
}
#endif


#if !CHANGE_BOUND_RIGHT
static Array evaluate(BoundBoth<Times> arg)
{
  return evaluateNumberBinary(arg, [](Number a, Number b) { return a * b; });
}
#endif


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


#if ADD_EACH
static Array evaluate(BoundRight<BoundOperator<First,Each>> /*arg*/)
{
  assert(false);
}
#endif


static Array evaluate(const char *arg)
{
  return makeCharArray(arg);
}


#if CHANGE_JOINING
static Array evaluate(Value value)
{
  return makeScalarArray(value);
}
#endif


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
  else if (arg1.shape.size() == 1 && arg2.shape.size() == 1) {
    Array a;
    a.shape.resize(1);
    a.shape[0] = 2;
    a.values.resize(2);
    a.values[0] = arg1;
    a.values[1] = arg2;
    return a;
  }
  else {
    cerr << "arg1.shape: " << arg1.shape << "\n";
    cerr << "arg2.shape: " << arg2.shape << "\n";
    assert(false);
  }
}


template <typename T>
static BoundRight<T> join(Primitive<T>, Array array)
{
  return { array };
}


#if CHANGE_JOINING
template <typename T>
static BoundRight<T> join(Primitive<T>, Value)
{
  assert(false);
}
#endif


#if CHANGE_JOINING
template <typename T>
static BoundRight<T> join(Primitive<T>, vector<Value> right)
{
  return { evaluate(right) };
}
#endif



#if !CHANGE_BOUND_RIGHT
template <typename T>
static BoundBoth<T> join(Array left, BoundRight<T> equal_array)
// rename equal_array
{
  return { left, equal_array.right };
}
#endif


#if CHANGE_JOINING
template <typename T>
static BoundBoth<T> join(Value left, BoundRight<T> equal_array)
{
  assert(false);
  // return { left, equal_array.right };
}
#endif


#if CHANGE_JOINING
template <typename T>
static BoundBoth<T> join(Value left, BoundBoth<T> equal_array)
{
  assert(false);
  // return { left, equal_array.right };
}
#endif


#if !CHANGE_BOUND_RIGHT
template <typename T>
static BoundBoth<T> join(Array arg1, BoundBoth<T> args)
{
  return { join(arg1, args.left), args.right };
}
#endif


static BoundRight<Equal> join(Primitive<Equal> left, BoundRight<Shape> right)
{
  return join(left, evaluate(right));
}


template <typename T>
static BoundRight<Operator<T>>
join(Operator<T> /*left*/, BoundRight<Iota> right)
{
  return {evaluate(right)};
}


template <typename T>
static Atop<Operator<T>,Primitive<Iota>>
join(Operator<T>, Primitive<Iota>)
{
  return {};
}


static Atop<BoundOperator<Plus,Reduce>,Primitive<Iota>>
join(Primitive<Plus>, Atop<Operator<Reduce>, Primitive<Iota>>)
{
  return {};
}


template <typename T>
static BoundRight<BoundOperator<Plus,T>>
join(Primitive<Plus> /*left*/, BoundRight<Operator<T>> right)
{
  return {right.right};
}


#if ADD_EACH
template <typename T>
static BoundRight<BoundOperator<First,T>>
join(Primitive<First> /*left*/, BoundRight<Operator<T>> right)
{
  return {right.right};
}
#endif


#if !CHANGE_BOUND_RIGHT
static BoundRight<Times>
join(Primitive<Times>, BoundBoth<Plus> right)
{
  return { evaluate(right) };
}
#endif


#if CHANGE_JOINING
static vector<Value> join(Value left, Value right)
{
  return { left, right };
}
#endif


#if CHANGE_EVALUATE_NUMBER
static Array join(Value left, Value right)
{
  Array a;
  a.shape.resize(1);
  a.shape[0] = 2;
  a.values.resize(2);
  a.values[0] = left;
  a.values[1] = right;
  return a;
}
#endif


#if CHANGE_JOINING
static vector<Value> join(Value left, vector<Value> right)
{
  vector<Value> result = { left };
  result.insert(result.end(), right.begin(), right.end());
  return result;
}
#endif


template <typename Arg>
static auto combine(Arg arg)
{
  return evaluate(arg);
}


template <typename Arg1, typename Arg2, typename ...Args>
static auto combine(Arg1 arg1, Arg2 arg2, Args ...args)
{
  return join(evaluate(arg1), combine(arg2, args...));
}


#if CHANGE_JOINING
static BoundRight<BoundOperator<Plus,Reduce>>
join(
  Atop<BoundOperator<Plus, Reduce>, Primitive<Iota> >,
  Value /*right*/
)
{
  assert(false);
  //return {evaluate(BoundRight<Iota>{right})};
}
#endif


static BoundRight<BoundOperator<Plus,Reduce>>
join(
  Atop<BoundOperator<Plus, Reduce>, Primitive<Iota> >,
  Array right
)
{
  return {evaluate(BoundRight<Iota>{right})};
}


namespace {
struct Placeholder {
  template <typename ...Args>
  auto operator()(Args ...args)
  {
    return evaluate(combine(args...));
  }

  static Primitive<Shape>  shape()  { return {}; }
  static Primitive<First>  first()  { return {}; }
  static Primitive<Equal>  equal()  { return {}; }
  static Primitive<Plus>   plus()   { return {}; }
  static Primitive<Times>  times()  { return {}; }
  static Primitive<Iota>   iota()   { return {}; }
  static Primitive<Reduce> reduce() { return {}; }
#if ADD_EACH
  static Primitive<Each>   each()   { return {}; }
#endif
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
  assert(_(_(_.plus, _.reduce, _.iota), 3) == _(6));
  assert(_(2, _.times, 3, _.plus, 1) == _(2*(3+1)));
  assert(_(_.first, 1,2,3) == _(1));
  assert(_(_.shape, _(1,2,3), _(4,5,6)) == _(2));
#if CHANGE_JOINING
  // We need to be able to differentiate between joining a scalar with
  // a list of scalars and joining a scalar with an array.
  assert(_(_.shape, 1, _(2,3)) == _(2));
#endif
#if ADD_EACH
  assert(_(_.first, _.each, 1, 2, 3, "ABC", _(9, 8, 7)) == _(1,2,3,'A',9));
#endif
}
