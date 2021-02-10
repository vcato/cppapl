#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include <memory>
#include <random>

#define ADD_OUTER 0

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

  friend ostream& operator<<(ostream &stream, Type type)
  {
    switch (type) {
      case Type::number:
        stream << "number";
        break;
      case Type::character:
        assert(false);
      case Type::array_ptr:
        stream << "array_ptr";
        break;
    }

    return stream;
  }

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

  explicit Value(const Value &arg)
  : type(arg.type)
  {
    createFrom(arg);
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
      cerr << "a.type: " << a.type << "\n";
      cerr << "b.type: " << b.type << "\n";
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
  explicit Array(const Array &) = default;
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
struct Keyword {
};
}


namespace {
template <typename T>
struct Function {
};
}


namespace {
template <typename T>
struct Operator {
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
#if ADD_OUTER
struct Outer;
#endif
struct Roll;
struct Replicate;
struct Empty;
struct Drop;
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


static Array makeArrayFromValues(Values v)
{
  if (v.empty()) {
    Array result;
    result.shape = {0};
    return result;
  }

  if (v.size() == 1) {
    return makeScalarArray(std::move(v[0]));
  }

  Array left;
  left.shape = { int(v.size()) };
  left.values = std::move(v);
  return left;
}


namespace {
struct Context {
  std::mt19937 random_engine;
};
}


static Value evaluate(int arg, Context &)
{
  return Number(arg);
}


static inline Value evaluate(Number arg, Context &)
{
  return arg;
}


static inline Value evaluate(char arg, Context &)
{
  return arg;
}


static Array evaluate(const char *arg, Context &)
{
  return makeCharArray(arg);
}


static Array evaluate(Value value, Context &)
{
  return makeScalarArray(std::move(value));
}


static Array evaluate(Values arg, Context &)
{
  return makeArrayFromValues(std::move(arg));
}


static Array evaluate(Atop<Function<Shape>,Array> arg, Context &)
{
  Array &array = arg.right;
  Array result;
  result.shape = { int(array.shape.size()) };

  for (int x : array.shape) {
    result.values.push_back(Number(x));
  }

  return result;
}


static Array evaluate(Atop<Function<Iota>, Array> arg, Context &)
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

  cerr << "arg.right: " << arg.right << "\n";
  assert(false);
}


static Array evaluate(Atop<Function<First>,Array> arg, Context &)
{
  if (arg.right.shape.size() == 0) {
    assert(false);
  }

  if (arg.right.shape.size() == 1) {
    return makeScalarArray(std::move(arg.right.values[0]));
  }

  assert(false);
}


static int roll(int value, Context &context)
{
  return std::uniform_int_distribution<int>(1,value)(context.random_engine);
}


namespace{
template <typename T>
struct Optional {
  bool has_value;

  Optional(T arg)
  : has_value(true),
    value(std::move(arg))
  {
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

  T &operator*()
  {
    assert(has_value);
    return value;
  }

  union {
    T value;
  };

  friend bool operator!=(const Optional &a, const Optional &b)
  {
    if (a.has_value != b.has_value) {
      assert(false);
    }

    if (a.has_value) {
      assert(b.has_value);
      return a.value != b.value;
    }

    assert(false);
  }
};
}


static Optional<int> maybeInteger(const Value &v)
{
  if (!v.isNumber()) {
    assert(false);
  }

  Number n = v.asNumber();
  int int_n = int(n);

  if (int_n != n) {
    assert(false);
  }

  return int_n;
}


static Array evaluate(Atop<Function<Roll>,Array> arg, Context &context)
{
  if (arg.right.shape.empty()) {
    if (arg.right.values[0].isNumber()) {
      int value = int(arg.right.values[0].asNumber());

      if (value != arg.right.values[0].asNumber()) {
        assert(false);
      }

      if (value > 0) {
        return makeScalarArray(Number(roll(value, context)));
      }
      assert(false);
    }
    assert(false);
  }

  if (arg.right.shape.size() == 1) {
    vector<Value> result;

    for (auto &x : arg.right.values) {
      if (!x.isNumber()) {
        assert(false);
      }

      int value = x.asNumber();

      if (value != arg.right.values[0].asNumber()) {
        assert(false);
      }

      result.push_back(Number(roll(value, context)));
    }

    return makeArrayFromValues(std::move(result));
  }

  cerr << "arg.right: " << arg.right << "\n";
  assert(false);
}


template <typename A, typename B, typename C>
static Atop<BoundOperator<A,B>,Function<C>>
evaluate(Atop<BoundOperator<A,B>,Function<C>> arg, Context &)
{
  return arg;
}


template <typename T>
static Function<T> evaluate(Function<T> (*)(), Context &)
{
  return {};
}


template <typename T>
static Operator<T> evaluate(Operator<T> (*)(), Context &)
{
  return {};
}


template <typename T>
static Keyword<T> evaluate(Keyword<T> (*)(), Context &)
{
  return {};
}


static Array evaluate(Keyword<Empty>, Context &)
{
  return makeArrayFromValues({});
}


static Array evaluate(Fork<Array,Function<Equal>,Array> arg, Context &)
{
  auto f = [](auto a, auto b){ return Number(a == b); };
  return evaluateBinary(std::move(arg.left), std::move(arg.right), f);
}


static Array
evaluate(Fork<Array,Function<Drop>,Array> arg, Context &/*context*/)
{
  if (arg.left.shape.size() == 0) {
    if (arg.right.shape.size() == 1) {
      Optional<int> maybe_n_to_drop = maybeInteger(arg.left.values[0]);

      if (!maybe_n_to_drop) {
        assert(false);
      }

      int n_to_drop = *maybe_n_to_drop;

      Array result;

      if (arg.right.shape[0] <= n_to_drop) {
        result.shape = {0};
        return result;
      }

      result.shape = {arg.right.shape[0] - n_to_drop};
      auto iter = arg.right.values.begin();
      iter += n_to_drop;

      for (; iter != arg.right.values.end(); ++iter) {
        result.values.push_back(std::move(*iter));
      }

      return result;
    }
  }

  assert(false);
}


template <typename T>
static Array evaluate(Fork<Values,T,Array> arg, Context &context)
{
  Array left = makeArrayFromValues(std::move(arg.left));
  using Return = Fork<Array, T, Array>;
  auto x = Return{std::move(left), arg.mid, std::move(arg.right)};
  return evaluate(std::move(x), context);
}


template <typename T, typename G>
static Array evaluateNumberBinary(Fork<Array,T,Array> arg, const G &g)
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

  return evaluateBinary(std::move(arg.left), std::move(arg.right), f);
}


static Array evaluate(Fork<Array,Function<Plus>,Array> arg, Context &)
{
  return
    evaluateNumberBinary(
      std::move(arg),
      [](Number a, Number b) { return a + b; }
    );
}


static Array evaluate(Fork<Array,Function<Times>,Array> arg, Context &)
{
  return
    evaluateNumberBinary(
      std::move(arg),
      [](Number a, Number b) { return a * b; }
    );
}


static Array evaluate(Array array, Context &)
{
  return array;
}


static Array evaluate(Atop<BoundOperator<Plus,Reduce>,Array> arg, Context &)
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


static Array evaluate(Atop<BoundOperator<First,Each>,Array> arg, Context &)
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

    return makeArrayFromValues(std::move(result));
  }
  else {
    assert(false);
  }

  assert(false);
}


template <typename T>
static Atop<Function<T>,Array> join(Function<T> left, Value right, Context &)
{
  return { left, makeScalarArray(std::move(right)) };
}


template <typename T>
static Atop<Function<T>,Array> join(Function<T> left, Array right, Context &)
{
  return { left, std::move(right) };
}


template <typename T>
static Atop<Function<T>,Array>
join(Function<T> left, Values right, Context &context)
{
  return { left, evaluate(std::move(right), context) };
}


template <typename T>
static Fork<Values,T,Array> join(Value left, Atop<T,Array> right, Context &)
{
  Values new_left;
  new_left.push_back(std::move(left));
  return { std::move(new_left), std::move(right.left), std::move(right.right) };
}


template <typename T>
static Fork<Values,T,Array>
join(Value left, Fork<Values,T,Array> right, Context &)
{
  Fork<Values,T,Array> result = std::move(right);
  result.left.insert(result.left.begin(), std::move(left));
  return result;
}


template <typename T1, typename T2>
static Atop<Function<T1>,Array>
join(Function<T1> left, Atop<Function<T2>,Array> right, Context &context)
{
  return join(left, evaluate(std::move(right), context), context);
}


template <typename T>
static Atop<Operator<T>,Array>
join(Operator<T> left, Atop<Function<Iota>,Array> right, Context &context)
{
  return { left, evaluate(std::move(right), context) };
}


static Atop<Atop<Operator<Product>,Function<Times>>, Array>
join(Operator<Product>, Atop<Function<Times>,Array> right, Context &)
{
  return {{},std::move(right.right)};
}


static Atop<Fork<Function<Plus>, Operator<Product>, Function<Times>>,Array>
join(
  Function<Plus>,
  Atop<Atop<Operator<Product>, Function<Times> >, Array> right,
  Context &
)
{
  return {{},std::move(right.right)};
}


#if ADD_OUTER
static Atop<Fork<Operator<Outer>, Operator<Product>, Function<Times>>,Array>
join(
  Operator<Outer> left,
  Atop<Atop<Operator<Product>, Function<Times> >, Array> right,
  Context&
)
{
  return {
    {
      std::move(left),
      std::move(right.left.left),
      std::move(right.left.right)
    },
    std::move(right.right)
  };
}
#endif


template <typename T>
static Atop<Operator<T>,Array>
join(Operator<T> left, Values right, Context &)
{
  return { left, makeArrayFromValues(std::move(right)) };
}


template <typename T>
static Atop<Operator<T>,Function<Iota>>
join(Operator<T>, Function<Iota>, Context &)
{
  return {};
}


static Atop<BoundOperator<Plus,Reduce>,Function<Iota>>
join(Function<Plus>, Atop<Operator<Reduce>, Function<Iota>>, Context &)
{
  return {};
}


template <typename T>
static Atop<BoundOperator<Plus,T>,Array>
join(Function<Plus> /*left*/, Atop<Operator<T>,Array> right, Context &)
{
  return { {}, std::move(right.right) };
}


template <typename T>
static Atop<BoundOperator<First,T>,Array>
join(Function<First> /*left*/, Atop<Operator<T>,Array> right, Context &)
{
  return {{},std::move(right.right)};
}


static Atop<Function<Times>,Array>
join(
  Function<Times>,
  Fork<Values,Function<Plus>,Array> right,
  Context &context
)
{
  return { {}, evaluate(std::move(right), context) };
}


static Array
evaluate(Fork<Array, Function<Replicate>, Array> arg, Context &)
{
  Array left = std::move(arg.left);
  Array right = std::move(arg.right);

  if (left.shape.empty() && right.shape.empty()) {
    Array result;

    if (!left.values[0].isNumber()) {
      assert(false);
    }

    Optional<int> maybe_n = maybeInteger(left.values[0]);

    if (!maybe_n) {
      assert(false);
    }

    int n = *maybe_n;

    for (int i=0; i!=n; ++i) {
      result.values.push_back(Value(right.values[0]));
    }

    result.shape = {int(result.values.size())};
    return result;
  }

  assert(false);
}


static Array
evaluate(
  Fork<
    Values,
    Fork<Function<Plus>,Operator<Product>,Function<Times>>,
    Array
  > arg
  , Context &
)
{
  Array left = makeArrayFromValues(std::move(arg.left));
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


static Values join(Value left, Value right, Context &)
{
  Values result;
  result.push_back(std::move(left));
  result.push_back(std::move(right));
  return result;
}


static Values join(Value left, Values right, Context &)
{
  Values result;
  result.push_back(std::move(left));

  for (auto &x : right) {
    result.push_back(std::move(x));
  }

  return result;
}


#if ADD_OUTER
template <typename T>
static Fork<Operator<Outer>,Operator<Product>,Function<T>>
join(Operator<Outer>, Atop<Operator<Product>,Function<T>>, Context &)
{
  assert(false);
}
#endif


template <typename Arg>
static auto combine(Context &context, Arg arg)
{
  return evaluate(std::move(arg), context);
}


template <typename Arg1, typename Arg2, typename ...Args>
static auto combine(Context &context, Arg1 arg1, Arg2 arg2, Args ...args)
{
  return
    join(
      evaluate(std::move(arg1), context),
      combine(context, std::move(arg2), std::move(args)...),
      context
    );
}


static Atop<BoundOperator<Plus,Reduce>,Array>
join(
  Atop<BoundOperator<Plus, Reduce>, Function<Iota> > left,
  Value right,
  Context &context
)
{
  return {
    left.left,
    evaluate(
      Atop<Function<Iota>,Array>
      {
        left.right, makeScalarArray(std::move(right))
      }
      , context
    )
  };
}


namespace {
struct Placeholder {
  Context context;

  template <typename ...Args>
  auto operator()(Args ...args)
  {
    return evaluate(combine(context, std::move(args)...), context);
  }

  static Function<Shape>     shape()     { return {}; }
  static Function<First>     first()     { return {}; }
  static Function<Equal>     equal()     { return {}; }
  static Function<Plus>      plus()      { return {}; }
  static Function<Times>     times()     { return {}; }
  static Function<Iota>      iota()      { return {}; }
  static Function<Roll>      roll()      { return {}; }
  static Function<Replicate> replicate() { return {}; }
  static Keyword<Empty>      empty()     { return {}; }
  static Function<Drop>      drop()      { return {}; }
  static Operator<Each>      each()      { return {}; }
  static Operator<Reduce>    reduce()    { return {}; }
  static Operator<Product>   product()   { return {}; }
#if ADD_OUTER
  static Operator<Outer>     outer()   { return {}; }
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
  assert(_(_.shape, _.shape, _.roll, 6) == _(0));
  assert(_(_.shape, _.roll, 6, 6) == _(2));
  assert(_(3, _.replicate, 4) == (_(4,4,4)));
  assert(_(1, _.drop, _.iota, 5) == (_(2,3,4,5)));
  assert(_(2, _.drop, _.iota, 5) == (_(3,4,5)));
  assert(_(2, _.drop, 4, 5) == _(_.empty));
  assert(_(3, _.drop, 4, 5) == _(_.empty));

#if ADD_OUTER
  assert(
    _(_(_.iota, 2), _.outer, _.product, _.times, _.iota, 2) ==
    _(2,2, _.reshape, 1, 2, 2, 4)
  );
#endif
}
