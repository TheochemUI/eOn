# Control-plane outcome. Trajectories stay readcon ConFrame files.
# JSON export is a decode of these three fields, not a second schema.
# Wire ordinals are API: never renumber existing fields.

@0xc3a91e70b4d25f18;

struct Outcome {
  jobType @0 :Text;
  statusCode @1 :Int32;
  statusText @2 :Text;
}
