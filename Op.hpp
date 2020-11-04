#pragma once

enum optype_t : uint32_t {
  DIAGONAL,
  OFF_DIAGONAL,
  IDENTITY,
};

struct Op {
  optype_t optype : 2;
  uint32_t opbond : 30;
};