When NEB `max_iterations` is not positive, `DynamicsSaddleSearch::run` passes the initial-band tangent to `MinModeSaddleSearch` instead of an empty mode.
