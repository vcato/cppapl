#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include <memory>

#define ADD_ROLL 0

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

  Value(Array arg);

  Value(Value &&arg)
  : type(arg.type)
  {
    createFrom(std::move(arg));
  }

  bool isArray() const { return type == Type::array_ptr; }
  bool isNumber() const { return type == Type::number; }
  bool isCharacter() const { return type == Type::character; }

  Number asNumber() const
  {
    assert(type == Type::number);
    return number;
  }

  char asCharacter() const
  {
    assert(type == Type::character);
    return character;
  }

  const Array &asArray() const
  {
    assert(type == Type::array_ptr);
    assert(array_ptr);
    return *array_ptr;
  }

  Array &asArray()
  {
    assert(type == Type::array_ptr);
    assert(array_ptr);
    return *array_ptr;
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
        return a.character == b.character;
      case Type::array_ptr:
        cerr << "a.array_ptr.get(): " << a.array_ptr.get() << "\n";
        cerr << "b.array_ptr.get(): " << b.array_ptr.get() << "\n";
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

  void createFrom(Value &&arg)
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
        new (&array_ptr) auto(std::move(arg.array_ptr));
        break;
    }
  }

  Value& operator=(Value arg)
  {
    destroy();
    type = arg.type;
    createFrom(std::move(arg));
    return *this;
  }

  ~Value() { destroy(); }
};
}


using Values = vector<Value>;


namespace {
struct Array {
  vector<int> shape;
  Values values;

  Array() = default;
  Array(const Array &) = delete;
  Array(Array &&) = default;

};
}


Value::Value(Array arg)
: type(Type::array_ptr), array_ptr(std::make_unique<Array>(std::move(arg)))
{
}


static ostream& operator<<(ostream& stream, const Value &v);
static ostream& operator<<(ostream& stream, const Array &a);


static ostream& operator<<(ostream& stream, const Value &v)
{
  if (v.isNumber()) {
    stream << v.asNumber();
  }
  else if (v.isCharacter()) {
    stream << "'" << v.asCharacter() << "'";
  }
  else if (v.isArray()) {
    stream << v.asArray();
  }
  else {
    assert(false);
  }

  return stream;
}


#if 1
static ostream& operator<<(ostream& stream, const Array &a)
{
  if (a.shape.empty()) {
    stream << a.values[0];
  }
  else if (a.shape.size() == 1) {
    stream << a.values;
  }
  else {
    stream << "shape: " << a.shape << "\n";
    assert(false);
  }

  return stream;
}
#endif


static bool operator==(const Array &a, const Array &b)
{
  if (a.shape != b.shape) {
    if (a.values.size() == 1 && b.values.size() == 1) {
      return a.values[0] == b.values[0];
    }

    cerr << "a: " << a << "\n";
    cerr << "b: " << b << "\n";
    assert(false);
  }

  if (a.values != b.values) {
    return false;
  }

  return true;
}
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
  result.values.push_back(std::move(value));
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
struct Function {
};
}


template <typename Left, typename Right>
struct Atop {
  Left left;
  Right right;
};


template <typename Left, typename Mid, typename Right>
struct Fork {
  Left left;
  Mid mid;
  Right right;
};


namespace {
struct Shape;
struct First;
struct Each;
struct Equal;
struct Plus;
struct Times;
struct Iota;
struct Reduce;
struct Product;
#if ADD_ROLL
struct Roll;
#endif
}


template <typename Function>
static Array evaluateBinary(Array left, Array right, Function f)
{
  if (left.values.size() == 1) {
    if (right.values.size() == 1) {
      Array result;
      result.shape = vector<int>{3};

      result.values.push_back(
        f(std::move(left.values[0]), std::move(right.values[0]))
      );

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
      *out++ = f(std::move(*in1++), std::move(*in2++));
    }

    return result;
  }

  cerr << "left.values: " << left.values << "\n";
  assert(false);
}


static Array makeArrayFromVector(Values v)
{
  if (v.empty()) {
    assert(false);
  }

  if (v.size() == 1) {
    return makeScalarArray(std::move(v[0]));
  }

  Array left;
  left.shape = { int(v.size()) };
  left.values = std::move(v);
  return left;
}


static Value evaluate(int arg)
{
  return Number(arg);
}


static inline Value evaluate(Number arg)
{
  return arg;
}


static inline Value evaluate(char arg)
{
  return arg;
}


static Array evaluate(Values arg)
{
  return makeArrayFromVector(std::move(arg));
}


static Array evaluate(Atop<Function<Shape>,Array> arg)
{
  Array &array = arg.right;
  Array result;
  result.shape = { int(array.shape.size()) };

  for (int x : array.shape) {
    result.values.push_back(Number(x));
  }

  return result;
}


static Array evaluate(Atop<Function<Iota>, Array> arg)
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


static Array evaluate(Atop<Function<First>,Array> arg)
{
  if (arg.right.shape.size() == 0) {
    assert(false);
  }

  if (arg.right.shape.size() == 1) {
    return makeScalarArray(std::move(arg.right.values[0]));
  }

  assert(false);
}


template <typename A, typename B, typename C>
static Atop<BoundOperator<A,B>,Function<C>>
evaluate(Atop<BoundOperator<A,B>,Function<C>> arg)
{
  return arg;
}


template <typename T>
static Function<T> evaluate(Function<T> (*)())
{
  return {};
}


template <typename T>
static Operator<T> evaluate(Operator<T> (*)())
{
  return {};
}


static Array evaluate(Fork<Values,Function<Equal>,Array> arg)
{
  auto f = [](auto a, auto b){ return Number(a == b); };
  Array left = makeArrayFromVector(std::move(arg.left));
  return evaluateBinary(std::move(left), std::move(arg.right), f);
}


template <typename T, typename G>
static Array evaluateNumberBinary(Fork<Values,T,Array> arg, const G &g)
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

  return
    evaluateBinary(
      makeArrayFromVector(std::move(arg.left)),
      std::move(arg.right),
      f
    );
}


static Array evaluate(Fork<Values,Function<Plus>,Array> arg)
{
  return
    evaluateNumberBinary(
      std::move(arg),
      [](Number a, Number b) { return a + b; }
    );
}


static Array evaluate(Fork<Values,Function<Times>,Array> arg)
{
  return
    evaluateNumberBinary(
      std::move(arg),
      [](Number a, Number b) { return a * b; }
    );
}


static Array evaluate(Array array)
{
  return array;
}


static Array evaluate(Atop<BoundOperator<Plus,Reduce>,Array> arg)
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


static Array evaluate(Atop<BoundOperator<First,Each>,Array> arg)
{
  if (arg.right.shape.empty()) {
    assert(false);
  }
  else if (arg.right.shape.size() == 1) {
    Values result;

    for (auto &x : arg.right.values) {
      if (x.isNumber()) {
        result.push_back(std::move(x));
      }
      else if (x.isArray()) {
        if (x.asArray().shape.empty()) {
          assert(false);
        }
        else if (x.asArray().shape.size() == 1) {
          if (x.asArray().values.empty()) {
            assert(false);
          }

          result.push_back(std::move(x.asArray().values[0]));
        }
        else {
          assert(false);
        }
      }
      else {
        assert(false);
      }
    }

    return makeArrayFromVector(std::move(result));
  }
  else {
    assert(false);
  }

  assert(false);
}


static Array evaluate(const char *arg)
{
  return makeCharArray(arg);
}


static Array evaluate(Value value)
{
  return makeScalarArray(std::move(value));
}


template <typename T>
static Atop<Function<T>,Array> join(Function<T> left, Value right)
{
  return { left, makeScalarArray(std::move(right)) };
}


template <typename T>
static Atop<Function<T>,Array> join(Function<T> left, Array right)
{
  return { left, std::move(right) };
}


template <typename T>
static Atop<Function<T>,Array> join(Function<T> left, Values right)
{
  return { left, evaluate(std::move(right)) };
}


template <typename T>
static Fork<Values,T,Array> join(Value left, Atop<T,Array> right)
{
  Values new_left;
  new_left.push_back(std::move(left));
  return { std::move(new_left), std::move(right.left), std::move(right.right) };
}


template <typename T>
static Fork<Values,T,Array> join(Value left, Fork<Values,T,Array> right)
{
  Fork<Values,T,Array> result = std::move(right);
  result.left.insert(result.left.begin(), std::move(left));
  return result;
}


static Atop<Function<Equal>,Array>
join(Function<Equal> left, Atop<Function<Shape>,Array> right)
{
  return join(left, evaluate(std::move(right)));
}


template <typename T>
static Atop<Operator<T>,Array>
join(Operator<T> left, Atop<Function<Iota>,Array> right)
{
  return { left, evaluate(std::move(right)) };
}


static Atop< Atop<Operator<Product>,Function<Times>>, Array>
join(Operator<Product>, Atop<Function<Times>,Array> right)
{
  return {{},std::move(right.right)};
}


static Atop<Fork<Function<Plus>, Operator<Product>, Function<Times>>,Array>
join(
  Function<Plus>,
  Atop<Atop<Operator<Product>, Function<Times> >, Array> right
)
{
  return {{},std::move(right.right)};
}


template <typename T>
static Atop<Operator<T>,Array>
join(Operator<T> left, Values right)
{
  return { left, makeArrayFromVector(std::move(right)) };
}


template <typename T>
static Atop<Operator<T>,Function<Iota>>
join(Operator<T>, Function<Iota>)
{
  return {};
}


static Atop<BoundOperator<Plus,Reduce>,Function<Iota>>
join(Function<Plus>, Atop<Operator<Reduce>, Function<Iota>>)
{
  return {};
}


template <typename T>
static Atop<BoundOperator<Plus,T>,Array>
join(Function<Plus> /*left*/, Atop<Operator<T>,Array> right)
{
  return { {}, std::move(right.right) };
}


template <typename T>
static Atop<BoundOperator<First,T>,Array>
join(Function<First> /*left*/, Atop<Operator<T>,Array> right)
{
  return {{},std::move(right.right)};
}


static Atop<Function<Times>,Array>
join(Function<Times>, Fork<Values,Function<Plus>,Array> right)
{
  return { {}, evaluate(std::move(right)) };
}


static Array
evaluate(
  Fork<
    Values,
    Fork<Function<Plus>,Operator<Product>,Function<Times>>,
    Array
  > arg
)
{
  Array left = makeArrayFromVector(std::move(arg.left));
  Array right = std::move(arg.right);

  if (left.shape.size() == 1 && right.shape.size() == 1) {
    if (left.values.size() == right.values.size()) {
      Number result = 0;
      size_t n = left.values.size();

      for (size_t i=0; i!=n; ++i) {
        assert(left.values[i].isNumber());
        assert(right.values[i].isNumber());
        result += left.values[i].asNumber()*right.values[i].asNumber();
      }

      return makeScalarArray(result);
    }
    assert(false);
  }
  assert(false);
}


static Values join(Value left, Value right)
{
  Values result;
  result.push_back(std::move(left));
  result.push_back(std::move(right));
  return result;
}


static Values join(Value left, Values right)
{
  Values result;
  result.push_back(std::move(left));

  for (auto &x : right) {
    result.push_back(std::move(x));
  }

  return result;
}


template <typename Arg>
static auto combine(Arg arg)
{
  return evaluate(std::move(arg));
}


template <typename Arg1, typename Arg2, typename ...Args>
static auto combine(Arg1 arg1, Arg2 arg2, Args ...args)
{
  return
    join(
      evaluate(std::move(arg1)),
      combine(std::move(arg2), std::move(args)...)
    );
}


static Atop<BoundOperator<Plus,Reduce>,Array>
join(
  Atop<BoundOperator<Plus, Reduce>, Function<Iota> > left,
  Value right
)
{
  return {
    left.left,
    evaluate(
      Atop<Function<Iota>,Array>
      {
        left.right, makeScalarArray(std::move(right))
      }
    )
  };
}


namespace {
struct Placeholder {
  template <typename ...Args>
  auto operator()(Args ...args)
  {
    return evaluate(combine(std::move(args)...));
  }

  static Function<Shape>   shape()   { return {}; }
  static Function<First>   first()   { return {}; }
  static Function<Equal>   equal()   { return {}; }
  static Function<Plus>    plus()    { return {}; }
  static Function<Times>   times()   { return {}; }
  static Function<Iota>    iota()    { return {}; }
#if ADD_ROLL
  static Function<Roll>    roll()    { return {}; }
#endif
  static Operator<Each>    each()    { return {}; }
  static Operator<Reduce>  reduce()  { return {}; }
  static Operator<Product> product() { return {}; }
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
  assert(_(3, _.equal, 3) == _(1));
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
  assert(_(_.shape, 1, _(2,3)) == _(2));
  assert(_(_.first, _.each, 1, 2, 3, "ABC", _(9, 8, 7)) == _(1,2,3,'A',9));
  assert(_(1,2,3,_.plus,_.product,_.times,4,5,6) == _(32));
#if ADD_ROLL
  assert(_(_.shape, _.shape, _.roll, 6) == _(0));
#endif
}
