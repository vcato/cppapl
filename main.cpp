#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <sstream>

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


template <typename T>
static void destroyObject(T &object)
{
  object.~T();
}


namespace {
struct Array;
}


using Values = vector<Array>;


namespace {
struct Box {
  vector<int> shape;
  vector<Array> values;

  friend bool operator==(const Box &a, const Box &b)
  {
    return a.shape == b.shape && a.values == b.values;
  }
};
}


namespace {
struct Boxed;
}


namespace {
class Array {
private:
  enum class Type { number, character, values , box };

  Type type;

  union {
    Number number;
    char character;
    Box box;
  };

  void createFrom(Array &&arg)
  {
    assert(type == arg.type);

    switch (type) {
      case Type::number:
        new (&number) auto(arg.number);
        break;
      case Type::character:
        new (&character) auto(arg.character);
        break;
      case Type::values:
      case Type::box:
        new (&box) auto(std::move(arg.box));
        break;
    }
  }

  void createFrom(const Array &arg)
  {
    assert(type == arg.type);

    switch (type) {
      case Type::number:
        new (&number) auto(arg.number);
        break;
      case Type::character:
        new (&character) auto(arg.character);
        break;
      case Type::values:
      case Type::box:
        new (&box) auto(arg.box);
        break;
    }
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
      case Type::values:
      case Type::box:
        destroyObject(box);
        break;
    }
  }

public:

  Array(Number arg)
  : type(Type::number), number(arg)
  {
  }

  Array(int arg)
  : type(Type::number), number(arg)
  {
  }

  Array(char arg)
  : type(Type::character), character(arg)
  {
  }

  Array(std::vector<int> shape, std::vector<Array> values)
  : type(Type::values),
    box{std::move(shape), std::move(values)}
  {
    assert(!this->shape().empty());
  }

  Array(Boxed arg);

  explicit Array(const Array &arg)
  : type(arg.type)
  {
    createFrom(arg);
  }

  Array(Array &&arg)
  : type(arg.type)
  {
    createFrom(std::move(arg));
  }

  ~Array() { destroy(); }

  const vector<int> &shape() const
  {
    assert(isNonScalar());
    return box.shape;
  }

  bool isNumber() const { return type == Type::number; }
  bool isCharacter() const { return type == Type::character; }
  bool isBox() const { return type == Type::box; }
  bool isNonScalar() const { return type == Type::values; }
  bool isSimple() const { return isNumber() || isCharacter(); }

  Number asNumber() const
  {
    assert(type == Type::number);
    return number;
  }

  char asCharacter() const
  {
    assert(false);
  }

  Values &asValues()
  {
    assert(isNonScalar());
    return box.values;
  }

  const Values &asValues() const
  {
    assert(isNonScalar());
    return box.values;
  }

  Array boxedValue() &&
  {
    assert(type == Type::box);
    assert(box.shape.empty());
    return std::move(box.values[0]);
  }

  Array& operator=(Array &&arg)
  {
    destroy();
    type = arg.type;
    createFrom(std::move(arg));
    return *this;
  }

  friend bool operator==(const Array &a, const Array &b)
  {
    if (a.type == Type::values && b.type == Type::number) {
      if (a.shape().size() == 1 && a.shape()[0] == 1) {
        return a.asValues()[0] == b;
      }

      assert(false);
    }

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
      case Type::values:
      case Type::box:
        return a.box == b.box;
    }
    assert(false);
  }

  friend ostream& operator<<(ostream &stream, Type type)
  {
    switch (type) {
      case Type::number:
        stream << "number";
        break;
      case Type::character:
        assert(false);
      case Type::values:
        stream << "non-scalar";
        break;
      case Type::box:
        assert(false);
        break;
    }

    return stream;
  }
};
}


namespace {
struct Boxed {
  Array array;
};
}


Array::Array(Boxed arg)
: type(Type::box),
  box{{}, {std::move(arg.array)}}
{
}


static ostream& operator<<(ostream& stream, const Array &a);


static void
printNonScalarArrayOn(
  ostream &stream, const vector<int> &shape, const Values &values
)
{
  if (shape.empty()) {
    stream << values[0];
  }
  else if (shape.size() == 1) {
    stream << values;
  }
  else if (shape.size() == 2) {
    stream << "[";

    for (int i=0; i!=shape[0]; ++i) {
      stream << " [";

      for (int j=0; j!=shape[1]; ++j) {
        stream << " " << values[i*shape[1] + j];
      }

      stream << " ]";
    }

    stream << " ]";
  }
  else {
    stream << "shape: " << shape << "\n";
    assert(false);
  }
}


static const Values &valuesOf(const Array &array)
{
  assert(array.isNonScalar());
  return array.asValues();
}


static Values &valuesOf(Array &array)
{
  assert(array.isNonScalar());
  return array.asValues();
}


static const vector<int> &shapeOf(const Array &a)
{
  return a.shape();
}


static ostream& operator<<(ostream& stream, const Array &v)
{
  if (v.isNumber()) {
    stream << v.asNumber();
  }
  else if (v.isCharacter()) {
    stream << "'" << v.asCharacter() << "'";
  }
  else if (v.isNonScalar()) {
    printNonScalarArrayOn(stream, shapeOf(v), valuesOf(v));
  }
  else {
    assert(false);
  }

  return stream;
}


namespace {
template <typename Function, typename Operator>
struct BoundOperator {
};
}


static Array makeCharArray(const char *arg)
{
  int n = strlen(arg);
  vector<int> shape = vector<int>{ n };
  Values values;

  for (int i=0; i!=n; ++i) {
    values.push_back(arg[i]);
  }

  assert(values[0].isCharacter());
  return Array(std::move(shape), std::move(values));
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
struct Context {
  std::mt19937 random_engine;
};
}


namespace {
template <typename F>
struct Expr {
  Context &context;
  F f;

  Expr(Context &context, F f)
  : context(context), f(std::move(f))
  {
  }

  bool evaluated = false;

  Expr(const Expr &) = delete;
  void operator=(const Expr &) = delete;
  void operator=(Expr &&) = delete;

  Expr(Expr &&arg)
  : context(arg.context),
    f(std::move(arg.f))
  {
    arg.evaluated = true;
  }

  ~Expr();
};
}


template <typename F>
static ostream& operator<<(ostream& stream, Expr<F> a)
{
  return stream << evaluateExpr(std::move(a));
}


namespace {
struct Var {
  Array *const ptr;
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
struct Reshape;
struct First;
struct Each;
struct Equal;
struct Plus;
struct Times;
struct Iota;
struct Reduce;
struct Product;
struct Outer;
struct Roll;
struct Replicate;
struct Enclose;
struct MemberOf;
struct Not;
struct Empty;
struct Assign;
struct Drop;
}


static bool isScalar(const Array &a)
{
  if (a.isNonScalar()) {
    assert(!a.shape().empty());
    return false;
  }
  else {
    return true;
  }
}


template <typename Function>
static Array evaluateBinary(Array left, Array right, Function f)
{
  if (isScalar(left)) {
    if (isScalar(right)) {
      return f(std::move(left), std::move(right));
    }

    assert(!right.shape().empty());
    Values values;

    for (auto &x : valuesOf(right)) {
      values.push_back(f(Array(left), std::move(x)));
    }

    return Array(right.shape(), std::move(values));
  }

  if (valuesOf(left).size() == 1) {
    if (valuesOf(right).size() == 1) {
      vector<int> shape = {1};
      Values values;

      values.push_back(
        f(std::move(valuesOf(left)[0]), std::move(valuesOf(right)[0]))
      );

      return Array(std::move(shape), std::move(values));
    }
    assert(false);
  }

  if (shapeOf(left) == shapeOf(right)) {
    vector<int> shape = left.shape();
    Values values;
    size_t n = valuesOf(left).size();
    auto in1 = valuesOf(left).begin();
    auto in2 = valuesOf(right).begin();

    for (size_t i=0; i!=n; ++i) {
      values.push_back(f(std::move(*in1++), std::move(*in2++)));
    }

    return Array(std::move(shape), std::move(values));
  }

  cerr << "left.values: " << valuesOf(left) << "\n";
  assert(false);
}


static Array makeArrayFromValues(Values v)
{
  if (v.empty()) {
    return Array({0}, {});
  }

  size_t n = v.size();

  if (n == 1) {
    return std::move(v[0]);
  }

  return Array({ int(n) }, std::move(v));
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

  explicit operator bool() const
  {
    return has_value;
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


static Optional<int> maybeInteger(const Array &a)
{
  if (isScalar(a)) {
    if (a.isNumber()) {
      Number n = a.asNumber();
      int int_n = int(n);

      if (int_n != n) {
        assert(false);
      }

      return int_n;
    }
    else {
      assert(false);
    }
  }
  else {
    return {};
  }
}


static int product(const vector<int> &arg)
{
  return std::accumulate(arg.begin(), arg.end(), 1, std::multiplies<>());
}


static void addNCopiesTo(Values &result_values, int n, const Array &v)
{
  assert(n >= 0);

  for (int i=0; i!=n; ++i) {
    result_values.push_back(Array(v));
  }
}


static void
replicateInto(
  Values &result_values,
  const Values &left_values,
  const Values &right_values,
  int n
)
{
  for (int i=0; i!=n; ++i) {
    Optional<int> maybe_count = maybeInteger(left_values[i]);

    if (!maybe_count) {
      assert(false);
    }

    addNCopiesTo(result_values, *maybe_count, right_values[i]);
  }
}


template <typename T, typename G>
static Array evaluateNumberBinary(Fork<Array,T,Array> arg, const G &g)
{
  auto f = [&](const Array& a, const Array& b)
  {
    if (a.isNumber() && b.isNumber()) {
      return g(a.asNumber(), b.asNumber());
    }
    else {
      assert(false);
    }
  };

  return evaluateBinary(std::move(arg.left), std::move(arg.right), f);
}


static Array evaluate(int arg, Context &)
{
  return Number(arg);
}


static Var evaluate(Var arg, Context &)
{
  assert(arg.ptr);
  return arg;
}


static inline Array evaluate(Number arg, Context &)
{
  return arg;
}


static inline Array evaluate(char arg, Context &)
{
  return arg;
}


static Array evaluate(Array value, Context &)
{
  return value;
}


static Array evaluate(Values arg, Context &)
{
  return makeArrayFromValues(std::move(arg));
}


static Array evaluate(Atop<Function<Shape>,Array> arg, Context &)
{
  Values values;

  if (isScalar(arg.right)) {
    return Array({0},{});
  }

  const vector<int> &shape = arg.right.shape();

  for (int x : shape) {
    values.push_back(Number(x));
  }

  return Array({ int(shape.size()) }, std::move(values));
}


static Array evaluate(Atop<Function<Iota>, Array> arg, Context &)
{
  if (isScalar(arg.right)) {
    Optional<int> maybe_n = maybeInteger(arg.right);

    if (!maybe_n) {
      assert(false);
    }

    int n = *maybe_n;
    assert(n >= 0);
    vector<int> shape = {n};

    Values values;

    for (int i=1; i<=n; ++i) {
      values.push_back(Number(i));
    }

    return Array(std::move(shape), std::move(values));
  }

  cerr << "arg.right: " << arg.right << "\n";
  assert(false);
}


static Array evaluate(Atop<Function<First>,Array> arg, Context &)
{
  if (isScalar(arg.right)) {
    if (arg.right.isBox()) {
      return std::move(arg.right).boxedValue();
    }
    assert(false);
  }

  if (arg.right.shape().size() == 1) {
    return std::move(valuesOf(arg.right)[0]);
  }

  assert(false);
}


static Array evaluate(Atop<Function<Roll>,Array> arg, Context &context)
{
  if (isScalar(arg.right)) {
    Optional<int> maybe_value = maybeInteger(arg.right);

    if (maybe_value) {
      int value = *maybe_value;

      if (value > 0) {
        return Number(roll(value, context));
      }

      assert(false);
    }

    assert(false);
  }

  if (arg.right.shape().size() == 1) {
    Values result;

    for (auto &x : valuesOf(arg.right)) {
      if (!x.isNumber()) {
        assert(false);
      }

      int value = x.asNumber();

      if (value != valuesOf(arg.right)[0].asNumber()) {
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
static Function<T> evaluate(Function<T> arg, Context &)
{
  return arg;
}


template <typename T>
static Operator<T> evaluate(Operator<T> arg, Context &)
{
  return arg;
}


template <typename T>
static Keyword<T> evaluate(Keyword<T> arg, Context &)
{
  return arg;
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
  if (isScalar(arg.left)) {
    if (arg.right.shape().size() == 1) {
      Optional<int> maybe_n_to_drop = maybeInteger(arg.left);

      if (!maybe_n_to_drop) {
        assert(false);
      }

      int n_to_drop = *maybe_n_to_drop;

      if (shapeOf(arg.right)[0] <= n_to_drop) {
        return Array({0}, {});
      }

      vector<int> shape = {shapeOf(arg.right)[0] - n_to_drop};
      auto iter = valuesOf(arg.right).begin();
      iter += n_to_drop;
      Values values;

      for (; iter != valuesOf(arg.right).end(); ++iter) {
        values.push_back(std::move(*iter));
      }

      return Array(std::move(shape), std::move(values));
    }
  }

  assert(false);
}


static Array evaluate(Fork<Array,Function<Reshape>,Array> arg, Context &)
{
  if (shapeOf(arg.left).size() != 1) {
    assert(false);
  }

  vector<int> shape;

  for (Array &x : valuesOf(arg.left)) {
    Optional<int> maybe_int_x = maybeInteger(x);

    if (!maybe_int_x) {
      assert(false);
    }

    shape.push_back(*maybe_int_x);
  }

  int n = product(shape);
  Values values;

  for (int i=0; i!=n; ++i) {
    values.push_back(valuesOf(arg.right)[i % valuesOf(arg.right).size()]);
  }

  return Array(std::move(shape), std::move(values));
}


static Array
evaluate(
  Fork<
    Array,
    Function<Fork<Operator<Outer>, Operator<Product>, Function<Times> > >,
    Array
  > arg,
  Context&
)
{
  if (arg.left.shape().size() != 1) {
    cerr << "arg.left: " << arg.left << "\n";
    cerr << "arg.left.shape.size(): " << arg.left.shape().size() << "\n";
    assert(false);
  }

  int n_rows = arg.left.shape()[0];

  if (arg.right.shape().size() != 1) {
    assert(false);
  }

  int n_cols = arg.right.shape()[0];

  vector<int> shape = { n_rows, n_cols };
  Values values;

  for (auto &x : valuesOf(arg.left)) {
    if (!x.isNumber()) {
      assert(false);
    }
  }

  for (auto &x : valuesOf(arg.right)) {
    if (!x.isNumber()) {
      assert(false);
    }
  }

  for (int i=0; i!=n_rows; ++i) {
    for (int j=0; j!=n_cols; ++j) {
      Number a = valuesOf(arg.left)[i].asNumber();
      Number b = valuesOf(arg.right)[j].asNumber();
      values.push_back(a * b);
    }
  }

  return Array(std::move(shape), std::move(values));
}


static bool isElementOf(const Array &a, const Array &b)
{
  assert(isScalar(a));

  if (isScalar(b)) {
    return (a == Array(b));
  }

  for (const Array &x : valuesOf(b)) {
    if (a == x) {
      return true;
    }
  }

  return false;
}


static Array evaluate(Fork<Array,Function<MemberOf>,Array> arg, Context &)
{
  if (isScalar(arg.left)) {
    Array left_value = std::move(arg.left);
    return isElementOf(left_value, arg.right);
  }

  size_t n = valuesOf(arg.left).size();
  Values values;

  for (size_t i=0; i!=n; ++i) {
    values.push_back(isElementOf(valuesOf(arg.left)[i], arg.right));
  }

  return Array(arg.left.shape(), std::move(values));
}


template <typename T>
static Array evaluate(Fork<Values,T,Array> arg, Context &context)
{
  Array left = makeArrayFromValues(std::move(arg.left));
  using Return = Fork<Array, T, Array>;
  auto x = Return{std::move(left), arg.mid, std::move(arg.right)};
  return evaluate(std::move(x), context);
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


static Array evaluate(Atop<BoundOperator<Plus,Reduce>,Array> arg, Context &)
{
  const Array &right = arg.right;

  if (right.shape().size() == 1) {
    Number total = 0;

    for (auto &x : valuesOf(right)) {
      if (!x.isNumber()) {
        assert(false);
      }

      total += x.asNumber();
    }

    return total;
  }

  assert(false);
}


static Array evaluate(Atop<BoundOperator<First,Each>,Array> arg, Context &)
{
  if (arg.right.shape().empty()) {
    assert(false);
  }
  else if (arg.right.shape().size() == 1) {
    Values result;

    for (auto &x : valuesOf(arg.right)) {
      if (x.isNumber()) {
        result.push_back(std::move(x));
      }
      else if (x.isNonScalar()) {
        if (shapeOf(x).empty()) {
          assert(false);
        }
        else if (shapeOf(x).size() == 1) {
          if (valuesOf(x).empty()) {
            assert(false);
          }

          result.push_back(std::move(valuesOf(x)[0]));
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


static Array evaluate(Fork<Array, Function<Replicate>, Array> arg, Context &)
{
  Array left = std::move(arg.left);
  Array right = std::move(arg.right);
  Values values;

  if (isScalar(left) && isScalar(right)) {
    Optional<int> maybe_n = maybeInteger(left);

    if (!maybe_n) {
      assert(false);
    }

    int n = *maybe_n;

    if (n < 0) {
      assert(false);
    }

    vector<int> shape = {n};
    Values values;

    for (int i=0; i!=n; ++i) {
      values.push_back(Array(right));
    }

    return Array(std::move(shape), std::move(values));
  }
  else if (shapeOf(left).size() == 1 && shapeOf(right).size() == 1) {
    if (shapeOf(left)[0] != shapeOf(right)[0]) {
      assert(false);
    }

    int n = shapeOf(left)[0];
    replicateInto(values, valuesOf(left), valuesOf(right), n);
  }
  else {
    assert(false);
  }

  size_t n = values.size();
  return Array({ int(n) }, std::move(values));
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

  if (left.shape().size() == 1 && right.shape().size() == 1) {
    if (valuesOf(left).size() == valuesOf(right).size()) {
      Number result = 0;
      size_t n = valuesOf(left).size();

      for (size_t i=0; i!=n; ++i) {
        assert(valuesOf(left)[i].isNumber());
        assert(valuesOf(right)[i].isNumber());
        result += valuesOf(left)[i].asNumber()*valuesOf(right)[i].asNumber();
      }

      return result;
    }
    assert(false);
  }
  assert(false);
}


static bool isNot(const Array &a)
{
  Optional<int> maybe_x = maybeInteger(a);

  if (maybe_x != 0 && maybe_x != 1) {
    assert(false);
  }

  return !*maybe_x;
}


static Array evaluate(Atop<Function<Not>, Array> arg, Context &)
{
  if (isScalar(arg.right)) {
    return isNot(arg.right);
  }

  size_t n = valuesOf(arg.right).size();
  Values values;

  for (size_t i=0; i!=n; ++i) {
    values.push_back(isNot(valuesOf(arg.right)[i]));
  }

  return Array(arg.right.shape(), std::move(values));
}


static Array evaluate(Atop<Function<Enclose>,Array> arg, Context &)
{
  if (arg.right.isSimple()) {
    return std::move(arg.right);
  }

  return Boxed{std::move(arg.right)};
}


template <typename T>
static Atop<Function<T>,Array> join(Function<T> left, Array right, Context &)
{
  return { left, std::move(right) };
}


static Atop<Keyword<Assign>,Array>
join(Keyword<Assign> left, Array right, Context &)
{
  return { left, std::move(right) };
}


template <typename T>
static auto join(T left, Var right, Context &context)
{
  assert(right.ptr);
  return join(left, Array(*right.ptr), context);
}


template <typename T>
static Atop<Function<T>,Array>
join(Function<T> left, Values right, Context &context)
{
  return { left, evaluate(std::move(right), context) };
}


static Array join(Var left, Atop<Keyword<Assign>, Array> right, Context&)
{
  assert(left.ptr);
  *left.ptr = std::move(right.right);
  return Array(*left.ptr);
}


static Atop<Keyword<Assign>,Array>
join(
  Keyword<Assign> left,
  Fork<Values, Function<Drop>, Array> right,
  Context& context
)
{
  return { std::move(left), evaluate(std::move(right), context) };
}


template <typename T>
static Fork<Values,T,Array> join(Array left, Atop<T,Array> right, Context &)
{
  Values new_left;
  new_left.push_back(std::move(left));
  return { std::move(new_left), std::move(right.left), std::move(right.right) };
}


template <typename T>
static Fork<Values,T,Array>
join(Array left, Fork<Values,T,Array> right, Context &)
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


static
Atop<
  Function<Fork<Operator<Outer>, Operator<Product>, Function<Times>>>,
  Array
>
join(
  Operator<Outer> /*left*/,
  Atop<Atop<Operator<Product>, Function<Times> >, Array> right,
  Context&
)
{
  return {
    {},
    std::move(right.right)
  };
}


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


template<typename T, typename U>
static Atop<Function<T>,Array>
join(
  Function<T> left,
  Fork<Values,Function<U>,Array> right,
  Context &context
)
{
  return { std::move(left), evaluate(std::move(right), context) };
}


static Values join(Array left, Array right, Context &)
{
  Values result;
  result.push_back(std::move(left));
  result.push_back(std::move(right));
  return result;
}


static Values join(Array left, Values right, Context &)
{
  Values result;
  result.push_back(std::move(left));

  for (auto &x : right) {
    result.push_back(std::move(x));
  }

  return result;
}


template <typename T>
static auto join(Var left, T right, Context &context)
{
  return join(Array(*left.ptr), std::move(right), context);
}


static Array combine(Context &, int arg)
{
  return Array(arg);
}


template <typename T>
static Function<T> combine(Context &, Function<T> arg)
{
  return arg;
}


template <typename T>
static Keyword<T> combine(Context &, Keyword<T> arg)
{
  return arg;
}


static Array combine(Context &, Array arg)
{
  return arg;
}


static Var combine(Context &, Var arg)
{
  return arg;
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


template <typename F>
static Expr<F> combine(Context&, Expr<F> expr)
{
  return expr;
}


template <typename T>
struct MakeRValue {
  T operator()(T &arg) const
  {
    return T(arg);
  }
};


template <typename F>
struct MakeRValue<Expr<F>> {
  Expr<F> operator()(Expr<F> &arg) const
  {
    return std::move(arg);
  }
};


template <typename...Args>
static auto evaluateInContext(Context &context, Args &&...args)
{
  return evaluate(combine(context, MakeRValue<Args>()(args)...), context);
}


template <typename F>
static auto evaluateExpr(Expr<F> &&expr)
{
  auto eval = [&](auto &&...args)
  {
    return
      evaluateInContext(
        expr.context,
        std::forward<decltype(args)>(args)...
      );
  };

  auto result = expr.f(eval);
  expr.evaluated = true;
  return result;
}


template <typename F>
Expr<F>::~Expr()
{
  if (!evaluated) {
    evaluateExpr(std::move(*this));
  }

  assert(evaluated);
}


template <typename F>
static auto evaluate(Expr<F> expr, Context &)
{
  return evaluateExpr(std::move(expr));
}


template <typename T, typename F>
static auto join(T left, Expr<F> right, Context &context)
{
  return join(std::move(left), evaluateExpr(std::move(right)), context);
}


static Atop<BoundOperator<Plus,Reduce>,Array>
join(
  Atop<BoundOperator<Plus, Reduce>, Function<Iota> > left,
  Array right,
  Context &context
)
{
  return {
    left.left,
    evaluate(
      Atop<Function<Iota>,Array>
      {
        left.right, std::move(right)
      }
      , context
    )
  };
}


template <typename T>
struct MakeRValue<const T&> {
  T operator()(const T &arg) const
  {
    return T(arg);
  }
};


template <typename T>
struct MakeRValue<T&> {
  Var operator()(T &arg) const
  {
    return Var{&arg};
  }
};


template <size_t n>
struct MakeRValue<const char (&)[n]> {
  Array operator()(const char *arg) const
  {
    return makeCharArray(arg);
  }
};


template <typename ...Args>
static auto makeExpr(Context &context, Args &&...args)
{
  auto lambda = [&](auto f){ return f(std::forward<decltype(args)>(args)...); };
  return Expr<decltype(lambda)>{context, lambda};
}


namespace {
struct Placeholder {
  Context context;

  template <typename ...Args>
  auto operator()(Args &&...args)
  {
    return makeExpr(context, std::forward<decltype(args)>(args)...);
  }

  static constexpr Function<Shape>     shape = {};
  static constexpr Function<Reshape>   reshape = {};
  static constexpr Function<First>     first = {};
  static constexpr Function<Equal>     equal = {};
  static constexpr Function<Plus>      plus = {};
  static constexpr Function<Times>     times = {};
  static constexpr Function<Iota>      iota = {};
  static constexpr Function<Roll>      roll = {};
  static constexpr Function<Replicate> replicate = {};
  static constexpr Function<MemberOf>  member_of = {};
  static constexpr Function<Not>       isnot = {};
  static constexpr Function<Drop>      drop = {};
  static constexpr Function<Enclose>   enclose = {};
  static constexpr Operator<Each>      each = {};
  static constexpr Operator<Reduce>    reduce = {};
  static constexpr Operator<Product>   product = {};
  static constexpr Operator<Outer>     outer = {};
  static constexpr Keyword<Empty>      empty = {};
  static constexpr Keyword<Assign>     assign = {};
};
}


constexpr Function<Shape>     Placeholder::shape;
constexpr Function<Reshape>   Placeholder::reshape;
constexpr Function<First>     Placeholder::first;
constexpr Function<Equal>     Placeholder::equal;
constexpr Function<Plus>      Placeholder::plus;
constexpr Function<Times>     Placeholder::times;
constexpr Function<Iota>      Placeholder::iota;
constexpr Function<Roll>      Placeholder::roll;
constexpr Function<Replicate> Placeholder::replicate;
constexpr Function<Drop>      Placeholder::drop;
constexpr Function<MemberOf>  Placeholder::member_of;
constexpr Function<Not>       Placeholder::isnot;
constexpr Function<Enclose>   Placeholder::enclose;
constexpr Keyword<Empty>      Placeholder::empty;
constexpr Keyword<Assign>     Placeholder::assign;
constexpr Operator<Each>      Placeholder::each;
constexpr Operator<Reduce>    Placeholder::reduce;
constexpr Operator<Product>   Placeholder::product;
constexpr Operator<Outer>     Placeholder::outer;


template <typename T>
static std::string str(const T &x)
{
  std::ostringstream stream;
  stream << x;
  return stream.str();
}


template <typename F>
static std::string str(Expr<F> expr)
{
  return str(evaluateExpr(std::move(expr)));
}


template <typename F>
static bool operator==(const Array &a, Expr<F> b)
{
  return a == evaluateExpr(std::move(b));
}


template <typename F1, typename F2>
static bool operator==(Expr<F1> a, Expr<F2> b)
{
  return evaluateExpr(std::move(a)) == evaluateExpr(std::move(b));
}


template <typename F>
static vector<int> shapeOf(Expr<F> expr)
{
  return shapeOf(evaluateExpr(std::move(expr)));
}


int main()
{
  Placeholder _;
  assert(_(1) == _(1));
  assert(!(_(2) == _(3)));
  assert(_(1,2) == _(1,2));
  assert(_(1,2,3) == _(1,2,3));
  assert(shapeOf(_(1,2,3)) == vector<int>{3});
  assert(_(_.shape, 1,2,3) == _(3));
  assert(_(3, _.equal, 3) == _(1));
  assert(_(3, _.equal, _.shape, 1,2,3) == _(1));
  assert(shapeOf(_(_.shape,1,2,3)) == vector<int>{1});
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
  assert(str(_(2,2, _.reshape, 1,2,3,4)) == "[ [ 1 2 ] [ 3 4 ] ]");

  assert(
    _(_(_.iota, 2), _.outer, _.product, _.times, _.iota, 2) ==
    _(2,2, _.reshape, 1, 2, 2, 4)
  );

  assert(_(1,0,2, _.replicate, 1,2,3) == _(1,3,3));

  {
    Array R = 5;
    _(R, _.assign, 1, _.drop, _.iota, 5);
    assert(R == _(2,3,4,5));
  }

  assert(_(1, _.member_of, 1) == _(1));
  assert(_(1, _.member_of, 1,2,3) == _(1));
  assert(_(1,2,3, _.member_of, 2) == _(0,1,0));
  assert(_(_.isnot, 1) == _(0));

  {
    Array R = 3;
    assert(_(_.iota, R) == _(1,2,3));
  }

  {
    Array R = 1;
    assert(_(R, _.plus, R) == _(2));
  }

  {
    Array R = 1;
    assert(_(R, _.plus, R, _.assign, 2) == _(4));
  }

  {
    Array R = 1;
    assert(_(_(R, _.plus, R), _.plus, R, _.assign, 2) == _(6));
  }

  {
    Array R = 5;
    assert(_(R, _.assign, 1, _.drop, _.iota, R) == _(2,3,4,5));
  }

  {
    Array R = 5;

    assert(
      _(
        _(_.isnot, R, _.member_of, R, _.outer, _.product, _.times, R),
        _.replicate, R, _.assign, 1, _.drop, _.iota, R
      ) == _(2,3,5)
    );
  }

  assert(_(_.shape, _.shape, _.enclose, 1, 2) == _(0));
  assert(_(_.first, _.enclose, 1,2) == _(1,2));
}
