---
layout: theme
---

# Blog Posts

{% for post in site.posts %}

#### {{ post.date | date: "%b %-d, %Y" }}

## [{{ post.title }}]({{ post.url }})

{% endfor %}

subscribe via [RSS]({{ "feed.xml" | prepend: site.baseurl }})

